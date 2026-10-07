# Figure 5: torque and orientation example

This example applies compatible Cartesian torques directly to an independently
integrated six-state spherical-eye model. It does not require participant data.

From the repository root in MATLAB:

```matlab
addpath('code/figure5');
reproduce_figure5;
```

The last PDF-margin step requires Python and `pypdf`. If Python is not on the
system path, pass its executable as `reproduce_figure5('/path/to/python')`.
MATLAB R2023b Update 4 was used to check this release. Results are written to
`outputs/figure5/`; rerunning does not modify `figures/manuscript/`.

| Setting | Value |
| --- | --- |
| Mass and radius | 0.0075 kg, 0.015 m |
| Spherical inertia | 6.75e-7 kg m² |
| Yaw/pitch proportional gains | 3725 / 4100 s⁻² |
| Derivative gains | 100 / 100 s⁻¹ |
| Reference time constant | 0.02 s |
| Final yaw/pitch | -28.7° / 10.8° |
| Initial state | Primary position, at rest |
| Integrator | ode45, RelTol=1e-10, AbsTol=1e-12, configured MaxStep=1e-3 s |
| Integration interval | 0–0.5 s; displayed interval 0–0.3 s |
| Saved sample spacing | 0.1 ms; this is not the adaptive internal integration step |

`feedback_states_corrected.csv` contains the six states in radians and radians/s.
`feedback_torques_corrected.csv` contains Cartesian torque in N m, ordered
torsion/x, pitch/y, yaw/z. `feedback_residuals_corrected.csv` gives the Listing
orientation, half-angle velocity, and torque-compatibility residuals.
`integration_checks.json` records tolerances, convergence and consistency checks.

Roll and roll rate are not projected onto Listing's relation during integration.
No torque deadband or roll-rate clamp is used. The maximum Listing orientation
residual in the checked run was approximately 3.45e-10 degrees. These checks
demonstrate numerical consistency of this example, not independent experimental
validation of human torsion or actuator performance.

The script retains an optional, disabled legacy force-allocation diagnostic.
It does not contribute to the applied torque or any result in this release.
