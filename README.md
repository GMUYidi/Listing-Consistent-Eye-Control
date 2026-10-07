# Listing-consistent eye dynamics and control: data and code

Research materials for **Listing-Consistent Reduced-Order Dynamics and Nonlinear
Feedback Control of Eye Movements**, by Yidi Huang, Qi Wei, Quentin Sanders, and
Ningshi Yao. This is the data/code release for the manuscript; it does not imply
journal acceptance. The manuscript has no assigned publication DOI.

Release **v1.0.0** contains the numerical eye-movement recordings underlying the
selected examples, retained analysis segments, fitting and simulation code,
archived manuscript figures, and documented rerun outputs.

## What is included

| Manuscript material | Released inputs and code | Reproduction status |
| --- | --- | --- |
| Figure 3: recorded trajectories | 9 yaw and 8 pitch segments; their 17 complete source records; `code/figures_3_4/` | Preserved observations and explicit display processing; new labels/layout |
| Figure 4: fitted poles | Same 17 segments and preserved fitting algorithm | Deterministic repository refit; the historical random state is unknown |
| Figure 5: torques and orientation | `code/figure5/`, parameters, states, torques, and numerical checks | Independently rerun six-state simulation matching the current manuscript example |
| Figure 6: six comparisons | Six complete source records, six retained windows, `code/figure6/` | Historical figure retained; new reconstruction explicitly documents unresolved historical alignment/window choices |
| Mathematical derivations | `code/theory/verify_math.py` | 14 symbolic identities and 200 synthetic consistency checks |

There are **23 distinct source recordings (36,754 saved rows)**, **17 fitting
segments (748 rows)**, and **six retained comparison windows (927 rows)**.
The full corresponding records are supplied so readers can inspect selection
and preprocessing. These are the records associated with this manuscript, not
the complete acquisition archive or a complete trial/exclusion log.

All fitting segments are from `participant_02`. Comparison panels (a), (c), (e),
and (f) are from `participant_01`; (b) and (d) are from `participant_02`.
Comparison records differ from fitting records, but selection was retrospective
and this is not a prospectively held-out or fully participant-independent test.
Participant and file identifiers are neutral labels.

## Start here

1. Download the [v1.0.0 release](https://github.com/GMUYidi/Listing-Consistent-Eye-Control/releases/tag/v1.0.0)
   or clone this repository and check out tag `v1.0.0`.
2. Run `python verify_release.py` to check file hashes, dataset counts, and the
   relationship between the full records and saved segments. This check needs
   only Python's standard library.
3. Install Python helpers with `python -m pip install -r requirements.txt`.
4. Run the components below, or inspect the committed results under `outputs/`.

From the repository root in MATLAB R2023b:

```matlab
addpath('code/figures_3_4');
reproduce_figures_3_4;

addpath('code/figure5');
reproduce_figure5;
% If needed: reproduce_figure5('/path/to/python')
```

From a terminal with Python 3.10 or newer:

```sh
python code/figure6/reconstruct_comparisons.py
python code/theory/verify_math.py
```

The MATLAB fitting driver uses Optimization Toolbox when available, with a
preserved `fminsearch` fallback. Signal Processing Toolbox also changes its
initial-guess routine if available. The checked release used both toolboxes;
alternative branches may produce different fitted poles. Figure 5 requires
MATLAB and Python/pypdf for PDF margins. Figure 6 uses NumPy, SciPy, and
Matplotlib. Exact checked versions and random seeds are recorded with the outputs.
No tracker hardware, SDK, private archive, or machine-specific path is required.

## Files and units

- `data/recordings/`: complete corresponding numerical recordings, with elapsed
  time, both SDK clocks, and eight raw/device-smoothed horizontal/vertical angle
  columns for the left and right eyes.
- `data/segments/`: retained fitting samples and comparison windows.
- `data/fitting_segments.csv` and `data/figure6_comparisons.csv`: record/segment
  correspondence, inclusive one-based source row indices, participant labels,
  processing information, and hashes.
- `figures/manuscript/`: unchanged figure assets currently used in the manuscript.
  Schematic Figures 1–2 have no experimental dataset and are not included here.
- `outputs/`: checked numerical outputs and rerun illustrations. Reruns overwrite
  this directory's component outputs, never the recorded data or archived figures.
- `metadata/checksums.json`: byte-level manifest of this release.
- `CITATION.cff`: citation metadata for this data/code release.

Recorded angles are in degrees; times are in seconds. Simulation CSV column
names specify radians, radians/s, or N m where applicable. Raw CSV angles retain
the SDK's horizontal/vertical convention. Manuscript coordinates use
**yaw = -horizontal angle** and **pitch = vertical angle**. Physical torque is
ordered `(torsion/x, pitch/y, yaw/z)`, while the independent controller inputs
are `(yaw/z, pitch/y)`.

“Raw” here refers to the saved SDK gaze-angle channel, not raw camera images.
Device-smoothed input values are preserved before additional analysis filtering.
All saved rows and duplicate SDK timestamps are retained. The 220 Hz software
polling target does not establish 220 independent image frames per second.
Elapsed time is shifted to the first saved row; both SDK clocks use one common
left-eye origin per recording, preserving their relative timing. See the
component READMEs for segment-specific time conventions.

## Scope and known limits

The release contains numerical observations without participant names, calendar
timestamps, eye videos, original archive filenames, or unrelated pupil/glint
metadata. No torsion measurement is supplied: the relevant saved torsion channels
contain no finite observations. Torsion curves derived from yaw and pitch are
Listing-law predictions, not independent torsion measurements.

The original stimulus geometry, calibration, device-filter settings, complete
trial/exclusion history, and historical optimization random state were not
recovered. The repository documents what is available instead of filling those
gaps by inference. In particular:

- Figure 4 reruns use seed `20261007`, which is a new documented repository seed.
  Their fitted poles need not be identical to the historical raster figure.
- The six Figure 6 record identities are supported by the retained source
  artwork, but the historical per-panel window, time/angle alignment, target
  selection, and complete batch renderer were not fully recovered. The provided
  reconstruction declares its conventions; it is not an exact recreation or a
  new validation of the archived figure. One retained window contains a long
  pre-movement interval, which remains visible without an inferred time shift.
- Figure 5 and symbolic checks test the equations and numerical implementation;
  they do not validate participant physiology or a physical actuator.

Detailed methods: [Figures 3–4](README_FIGURES_3_4.md),
[Figure 5](README_FIGURE_5.md), [Figure 6](README_FIGURE_6.md).
Use the versioned release and `CITATION.cff` when referring to these materials.
