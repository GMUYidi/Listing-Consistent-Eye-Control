# Data and analysis supporting Figures 3 and 4

This package contains the 17 selected left-eye segments used for response fitting
(9 yaw and 8 pitch), the 17 corresponding complete saved numerical recordings,
and the fitting and plotting code. All selected fitting segments belong to
`participant_02`. The participant labels are neutral identifiers, not names.
The selection was retrospective. These are not the complete candidate population
or an exclusion log, and the number of source files is not a count of all trials.

## Run

In MATLAB, add `code/figures_3_4` to the path, then call:

```matlab
reproduce_figures_3_4
```

The driver locates the data relative to itself. It writes only to
`outputs/figures_3_4` by default; an alternative output directory may be passed
as its argument. It produces CSV fitted parameters, CSV display curves,
PDF/PNG replots, and a JSON environment/seed record. Figures are created invisibly.
To redraw existing results without refitting, call
`reproduce_figures_3_4(output_directory, true)` after a completed run.

The code targets MATLAB R2023b. The preserved fitter uses `lsqcurvefit` when
Optimization Toolbox is available and a projected `fminsearch` fallback otherwise.
Its starting-value estimator uses `findpeaks` when available (Signal Processing
Toolbox); otherwise it uses the retained maximum-based fallback. Toolbox choices
can change numerical results. No ViewPoint SDK, device, Python, or original archive
is needed to run the public MATLAB analysis.

## Historical figures versus repository reruns

The driver uses the explicitly documented repository seed `20261007` with MATLAB's
`twister` generator. This is a new deterministic rerun seed; the random seed and
exact optimization environment that produced the manuscript's historical pole
plots were not recorded. The generated Figure 4 replots are labeled repository
refits and are not claimed to recover every historical pole or image pixel.
Running this package does not overwrite or revise manuscript figures or gains.

The Figure 3 display processing follows the retained plotting scripts, with
neutral labels and a new layout. It is not the input resampling used for fitting:
fitting retains the measured, nonuniform saved times. Figure 3 display curves
are explicitly derived visualization samples and are not extra observations.

## Files and time conventions

- `data/fitting_segments.csv` identifies each segment, its complete-recording CSV,
  participant, fitted axis, inclusive one-based source row indices, preprocessing
  windows, measured amplitude, timestamp counts, and public-file hashes.
- `data/segments/yaw_01.csv` through `yaw_09.csv` and `pitch_01.csv` through
  `pitch_08.csv` contain `sample_index`, `time_s`, `angle_deg`, and
  `left_sdk_time_s`. `time_s` is elapsed time relative to the original manually
  selected window origin; its first value need not be zero. `angle_deg` is the
  saved device-smoothed left-eye channel after the sign conversion described below,
  before the fitting routine's additional moving average. The segment SDK time
  is relative to that segment's first left-eye SDK timestamp and is not used in
  the fitter.
- `data/recordings/fit_yaw_01.csv` etc. retain the complete corresponding numerical
  record, including both eyes and the raw and device-smoothed SDK gaze-angle
  channels. Their `elapsed_time_s` starts at the first saved elapsed-time value.
  Both SDK columns have the **same** origin: the first finite left-eye SDK timestamp
  in that recording. Thus the left/right SDK timing relation is retained. The
  elapsed clock and SDK clock are distinct, and their zero points must not be
  interpreted as synchronized image acquisition times.

To reconstruct a segment from its complete recording, select the inclusive rows
`source_start_sample_1based:source_end_sample_1based`, then subtract
`segment_time_origin_record_s` from the recording's elapsed time. The selected
window's requested end is `selected_window_end_record_s`. The angles are in
degrees. The source columns retain the SDK horizontal/vertical semantics:
**manuscript yaw = negative horizontal angle; pitch = vertical angle**.
For the retained source columns this is, for example,
`yaw = -smoothed_horizontal_left_deg` and
`pitch = smoothed_vertical_left_deg`. Raw and smoothed values are retained
without changing their signs in the complete-recording CSVs.

Here “raw” means the SDK `GazeAngle` numeric channel, not eye images or a guarantee
of entirely unprocessed sensor measurements. The recorded device smoothing
settings were not recovered. CSVs contain no names, archive filenames, original
MAT headers, calendar times, video, gaze images, or unrelated metadata. No
torsion measurement is supplied because these records have no finite saved
torsion observations. No gaze values have been imputed in this export.

The 17 selected windows contain 748 saved polling rows and 483 distinct left-eye
SDK timestamps counted within windows; 265 adjacent timestamp repeats are retained.
A 220 Hz software polling target is not proof of 220 independent images per second.
The complete records permit checking window selection and numerical preprocessing;
they do not supply missing calibration, stimulus geometry, or an exclusion history.

## Preserved fitting procedure

The numerical core is adapted from the retained fitter. Only file loading,
private defaults, commented plotting, console reporting, and a workspace side
effect were removed or replaced. The response fitting calculations are preserved.

1. Apply `movmean` with `max(3,round(0.02*n))` samples when `n>15`, using shortened
   boundary windows. This equals three samples for every selected segment.
2. Use the first filtered sample as baseline and the median from one-based sample
   `round(0.85*n)` through `n` as the endpoint. The initial 15% median is used only
   to determine the response direction.
3. Fit from three samples before the first 5% displacement crossing (bounded by
   the segment start) through the first 98% crossing, inclusive.
4. The selected output branch is underdamped. Stage 1 fixes baseline and fits
   damping ratio, natural frequency, onset shift, and amplitude. Stage 2 fixes
   the measured baseline and amplitude and refines the other three parameters.
   Each stage retains the best of 20 starts, without a guarantee of global optimality.
5. Bounds are `0.02 <= zeta <= 0.98`, `0.001 <= omega_n <= 10000` per second, and
   `t_first - 0.2*(t_last-t_first) <= t0 <= t_last`, using the **full saved segment**
   for these time bounds. Stage 1 amplitude is unbounded in the retained routine;
   finite random-start bounds use +/-1000 for otherwise infinite bounds. The
   solver limits are 50,000 function evaluations and 2,000 iterations per start.

The retained function also computes two overdamped stages before the underdamped
stages, even when underdamped output is explicitly selected. This sequence is
preserved because it affects random-number consumption. It does **not** mean the
manuscript selected underdamping by comparison with the overdamped model.
The output RMSE and MAE are evaluated over all filtered segment samples; the
optimization objective uses only the fitting window. Gains are `kp=omega_n^2`
and `kd=2*zeta*omega_n`. The manuscript's representative engineering gains are
not replaced by an automatic average of these rerun parameters.

## Display-only procedure for Figure 3

The retained display code subtracts the first saved time and angle, crops to
180 ms for yaw or 140 ms for pitch, interpolates within the saved range with PCHIP
on a 1 ms grid, applies MATLAB Gaussian smoothing with a maximum odd window of
29 samples, and subtracts the first smoothed value again. The final displayed
value is held constant to the display endpoint if needed. The retained initial
offset rules are 5 degrees for yaw and 10 degrees for pitch. The driver checks
these rules for the selected set and fails rather than silently dropping a segment.
These display operations neither determine the fit's amplitude nor create
additional samples for fitting. Neutral segment labels replace archive labels;
they do not denote stimulus amplitudes.
