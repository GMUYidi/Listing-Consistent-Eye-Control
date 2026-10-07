"""Reconstruct comparisons from released recordings; not the historical Fig. 6 renderer.

This evaluates the ideal reduced closed-loop equations used by the manuscript.
Full six-state torque integration is provided separately with Figure 5.
No model parameter, sample offset, or time shift is optimized here.
"""
from pathlib import Path
import argparse
import csv
import json
import platform
import numpy as np
import scipy
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
KP = np.array([3725.0, 4100.0])
KD = np.array([100.0, 100.0])
TAU_REF = 0.02


def read_numeric(path):
    # Python's float reader preserves the exported double-precision values.
    with path.open(encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    return {key: np.array([float(row[key]) for row in rows]) for key in rows[0]}


def listing_proxy_deg(yaw_deg, pitch_deg):
    yaw, pitch = np.deg2rad(yaw_deg), np.deg2rad(pitch_deg)
    if np.any(np.abs(yaw) >= np.pi / 2) or np.any(np.abs(pitch) >= np.pi / 2):
        raise ValueError("The chosen coordinates leave the stated Listing/Euler domain.")
    return np.rad2deg(2 * np.arctan(np.tan(yaw / 2) * np.tan(pitch / 2)))


def validate_segment(record, segment, meta):
    ix = segment["source_row_1based"].astype(int) - 1
    if not np.array_equal(ix + 1, segment["source_row_1based"]):
        raise ValueError("Nonintegral source index")
    if ix[0] != int(meta["source_start_row_1based"]) - 1 or ix[-1] != int(meta["source_end_row_1based"]) - 1:
        raise ValueError("Source indices disagree with manifest")
    for key in record:
        if not np.array_equal(record[key][ix], segment[key], equal_nan=True):
            raise ValueError(f"Segment differs from source recording: {key}")
    expected = segment["elapsed_time_s"] - float(meta["segment_time_origin_record_s"])
    if not np.allclose(expected, segment["segment_time_s"], atol=2e-14, rtol=0):
        raise ValueError("Segment time origin disagrees with manifest")


def integrate(target_deg, end_time):
    target = np.deg2rad(target_deg)

    def rhs(t, state):
        q, qdot = state[[0, 2]], state[[1, 3]]
        # Legacy filtered-position PD convention: no reference-derivative feedforward.
        reference = target * (1 - np.exp(-t / TAU_REF))
        acceleration = -KP * (q - reference) - KD * qdot
        return [qdot[0], acceleration[0], qdot[1], acceleration[1]]

    result = solve_ivp(rhs, (0.0, end_time), np.zeros(4), method="RK45",
                       rtol=1e-10, atol=1e-12, max_step=1e-3, dense_output=True)
    if not result.success:
        raise RuntimeError(result.message)
    return result


def write_numeric(path, columns):
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        writer.writerows(zip(*columns.values()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs/figure6_reconstruction")
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    with (ROOT / "data/figure6_comparisons.csv").open(encoding="utf-8", newline="") as f:
        manifest = list(csv.DictReader(f))
    if len(manifest) != 6:
        raise ValueError("Expected six comparison records")
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    report = []
    for ax, meta in zip(axes.flat, manifest):
        record = read_numeric(ROOT / meta["source_record_csv"])
        segment = read_numeric(ROOT / meta["selected_segment_csv"])
        validate_segment(record, segment, meta)
        time = segment["segment_time_s"]
        # These are componentwise gaze changes, not a recovered 3-D frame rotation.
        yaw = -segment["smoothed_horizontal_left_deg"]
        pitch = segment["smoothed_vertical_left_deg"]
        yaw, pitch = yaw - yaw[0], pitch - pitch[0]
        if not np.isfinite(np.c_[time, yaw, pitch]).all():
            raise ValueError("Missing selected left-eye data; no implicit interpolation is allowed")
        tail = int(meta["reconstruction_tail_samples"])
        target = np.array([np.median(yaw[-tail:]), np.median(pitch[-tail:])])
        recorded_target = np.array([float(meta["reconstruction_yaw_target_deg"]),
                                    float(meta["reconstruction_pitch_target_deg"])])
        if not np.allclose(target, recorded_target, atol=1e-12, rtol=0):
            raise ValueError("Target extraction differs from manifest")
        solution = integrate(target, float(time[-1]))
        state = solution.sol(time)
        sim_yaw, sim_pitch = np.rad2deg(state[[0, 2]])
        sim_roll = listing_proxy_deg(sim_yaw, sim_pitch)
        measured_proxy = listing_proxy_deg(yaw, pitch)
        write_numeric(out / f"{meta['comparison_id']}_reconstructed.csv", {
            "segment_time_s": time, "source_row_1based": segment["source_row_1based"].astype(int),
            "measured_yaw_change_deg": yaw, "measured_pitch_change_deg": pitch,
            "listing_proxy_from_centered_gaze_deg": measured_proxy,
            "simulated_yaw_deg": sim_yaw, "simulated_pitch_deg": sim_pitch,
            "simulated_listing_roll_deg": sim_roll,
        })
        fine_time = np.linspace(0.0, float(time[-1]), 1001)
        fine_state = solution.sol(fine_time)
        sy, sp = np.rad2deg(fine_state[[0, 2]])
        sr = listing_proxy_deg(sy, sp)
        ax.plot(time, yaw, color="#bd4740", lw=1.6, label="Measured yaw change")
        ax.plot(time, pitch, color="#246d9e", lw=1.6, label="Measured pitch change")
        ax.plot(time, measured_proxy, color="#7c7c7c", ls=":", lw=1.4, label="Listing proxy from centered gaze")
        ax.plot(fine_time, sy, color="#bd4740", ls="--", lw=1.6, label="Reconstructed model yaw")
        ax.plot(fine_time, sp, color="#246d9e", ls="--", lw=1.6, label="Reconstructed model pitch")
        ax.plot(fine_time, sr, color="#20884c", ls="--", lw=1.6, label="Model Listing roll")
        ax.set(title=f"({meta['figure6_panel']}) {meta['comparison_id']} / {meta['participant_id']}",
               xlabel="Saved segment time (s)", ylabel="Angle / component change (deg)")
        ax.grid(True, alpha=0.25)
        ax.margins(x=0.01, y=0.1)
        report.append({"comparison_id": meta["comparison_id"], "participant_id": meta["participant_id"],
            "source_samples": len(record["source_row_1based"]), "selected_samples": len(time),
            "source_segment_numeric_equality": True, "target_yaw_deg": float(target[0]),
            "target_pitch_deg": float(target[1]), "tail_samples": tail,
            "solver_function_evaluations": solution.nfev,
            "actual_max_accepted_step_s": float(np.max(np.diff(solution.t))),
            "applied_extra_time_shift_s": 0.0, "historical_figure6_exact_reproduction": False})
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncols=3, fontsize=9)
    fig.suptitle("Reconstruction from retained segments (historical Figure 6 processing is not fully recovered)", fontsize=14)
    fig.savefig(out / "comparison_reconstruction.png", dpi=180)
    plt.close(fig)
    summary = {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__, "solver": "scipy.integrate.solve_ivp/RK45",
        "rtol": 1e-10, "atol": 1e-12, "max_step_s": 1e-3,
        "kp_s_minus2": KP.tolist(), "kd_s_minus1": KD.tolist(), "tau_ref_s": TAU_REF,
        "initial_yaw_pitch_deg": [0.0, 0.0], "initial_yaw_pitch_rates_deg_s": [0.0, 0.0],
        "target_rule": "median of final ceil(15% of saved samples), at least 3, minus first saved value",
        "comparison_records": report}
    (out / "reconstruction_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
