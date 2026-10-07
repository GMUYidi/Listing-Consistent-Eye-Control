# Figure 6 data and transparent reconstruction

This release contains six device-output recordings associated with panels (a)-(f) of the article's Figure 6, and one existing saved segment per recording. Their panel identities are supported by the retained source artwork. Four records are from `participant_01` and two from `participant_02`. The separate 17-record fitting set is entirely from `participant_02`; the six comparisons use distinct recordings, but do not constitute a prospectively held-out or wholly participant-independent validation set.

The historical plot's full processing recipe was not recovered. The new script therefore produces an explicitly labeled reconstruction. It does not overwrite or purport to regenerate the archived article figure exactly.

## Files and numerical provenance

- `data/recordings/comparison_01.csv` through `comparison_06.csv`: all rows of the six source recordings, restricted to the timestamps and eight raw/device-smoothed angle columns used for reproducibility.
- `data/segments/comparison_01.csv` through `comparison_06.csv`: exact rows from existing saved segment files, with both eyes and their original segment time retained.
- `data/figure6_comparisons.csv`: participant/panel mapping, inclusive source row indices, saved window bounds, sample counts, reconstruction target rule, and SHA-256 hashes of the public CSV files.
- `code/figure6/reconstruct_comparisons.py`: portable reconstruction and source-to-segment verification.

No new manual windows were selected. When multiple saved variants existed, a retained longer combined yaw/pitch segment was used. For panels (b), (c), and (e), these longer segments were retained with the model files; for (a), (d), and (f), they were retained with the acquisition data. This choice makes the available trajectory transparent; it does not establish which variant was displayed in the historical figure. Row counts for (a)-(f) are 110, 123, 110, 127, 126, and 331.

CSV floating-point values use round-trip decimal serialization. All retained angle values, sample ordering, and duplicate SDK times are preserved. Only time origins and field/record identifiers are changed. The source MAT's unrelated pupil, glint, and screen gaze-point channels are omitted. Its torsion channels contain no finite measurements for these records, so no independent measured torsion is released or claimed.

## CSV columns

`source_row_1based` is the 1-based row in the full source recording; it is not reset in the selected-segment CSV. `elapsed_time_s` equals acquisition elapsed time minus its first recorded value. `sdk_time_left_s` and `sdk_time_right_s` both subtract the same value: the first finite left-eye SDK timestamp in that recording. This preserves repeated timestamps and the timing relation between eyes. The CSVs retain every acquisition row; a nominal 220 Hz acquisition loop does not mean that every row is an independent device frame.

The following eight columns preserve the original device horizontal/vertical angles in degrees:

```
raw_horizontal_left_deg       raw_vertical_left_deg
raw_horizontal_right_deg      raw_vertical_right_deg
smoothed_horizontal_left_deg  smoothed_vertical_left_deg
smoothed_horizontal_right_deg smoothed_vertical_right_deg
```

Here `raw` means the SDK's `GazeAngle_cell` output, not camera images or unprocessed optical observations. `smoothed` means `GazeAngleSmoothed_cell`; the acquisition-specific device smoothing settings are not recovered. Manuscript coordinates use `yaw = -horizontal` and `pitch = vertical`. The release does not add smoothing, deduplicate, interpolate, or resample the measured source columns.

Selected segments additionally include `segment_time_s`, copied from the saved segment's `time_pad`. The manifest's `segment_time_origin_record_s` is the saved manual start bound relative to the first record row, so `segment_time_s = elapsed_time_s - segment_time_origin_record_s` up to floating-point subtraction roundoff. The bounds describe an archived manual selection; they are not established physiological onset/offset times or programmed stimulus directions.

## Run the reconstruction

From the repository root, using Python 3.10 or newer:

```sh
python -m pip install -r code/figure6/requirements.txt
python code/figure6/reconstruct_comparisons.py
```

Paths are resolved from the script location, so the second command can also be invoked from another working directory. It writes one six-panel PNG, six derived CSVs, and a run summary to `outputs/figure6_reconstruction/`. Use `--output-dir PATH` to choose another output directory. Tested versions are recorded in the run summary. Measured source CSVs are never changed.

The driver verifies that each segment equals the indicated rows of its full source recording. It uses left-eye device-smoothed angles, converts horizontal to yaw, then subtracts each component's first saved value. This componentwise centering is a transparent plotting convention, not a recovered three-dimensional calibration or reference-frame rotation. No extra time shift is fitted; saved `segment_time_s` is used directly, and model integration begins at zero with zero yaw/pitch and rates.

For the new reconstruction, each target is the median of the final 15% of saved samples (rounded up, minimum three) minus the first value. This explicit rule is not asserted to be the historical target-selection rule. The ideal reduced dynamics are

```
reference(t) = target * (1 - exp(-t / 0.02))
q_ddot = -KP * (q - reference(t)) - KD * q_dot
KP = diag(3725, 4100) s^-2; KD = diag(100, 100) s^-1
```

They use the filtered-position PD convention of the numerical example, without reference-derivative feedforward. Integration uses SciPy RK45 with relative tolerance 1e-10, absolute tolerance 1e-12, and maximum step 0.001 s. This reproduces the ideal controlled yaw/pitch equations, not a separate muscle-force or full six-state numerical validation; the latter is provided with Figure 5.

Model roll and a separately labeled proxy from centered gaze are evaluated using `2*atan(tan(yaw/2)*tan(pitch/2))`. The proxy is computed from the centered gaze components and the assumed Listing relation; it is not measured torsion. The plotted model roll is likewise reconstructed algebraically. The selected observations are drawn at their original irregular timestamps; dense evaluation is used only for the simulated curves.

## Historical reproduction limits

The recovered files establish source-record and saved-window provenance. They do not fully establish the historical figure's per-panel chosen window, additional time/angle offsets, interpolation or display smoothing, target-selection rule, or complete per-panel simulation configuration. No complete original six-panel batch renderer was recovered. The new script makes these choices explicit and does not fit them to the archived image. Original stimulus directions, calibration details, screen geometry, and device filter settings cannot be inferred from these CSVs.

In particular, the saved `comparison_06` window includes approximately 0.4 s before its large movement. The model starts at saved segment time zero under the declared convention. The resulting time mismatch is visible in the reconstruction; it is not removed by silently selecting a different window or estimating an onset shift.
