function report = verify_dynamic_model_integration
% Run the plotting script and restore the user's figure settings afterwards.
previous_visibility = get(groot, 'defaultFigureVisible');
existing_figures = findall(groot, 'Type', 'figure');
cleanup = onCleanup(@() restore_figures(previous_visibility, existing_figures)); %#ok<NASGU>
set(groot, 'defaultFigureVisible', 'off');
report = run_integration_checks();
end

function report = run_integration_checks
% The original script's clear command affects this workspace only.
run(fullfile(fileparts(mfilename('fullpath')), 'Dynamic_model1210.m'));
verify_dir = output_dir;
saved_result = load(fullfile(verify_dir,'feedback_result_corrected.mat'),'params');
params = saved_result.params; % Explicit binding avoids a params.m path collision.
assert(m == 0.0075 && r == 0.015);
assert(abs(t(end)-T_max) < eps && numel(t) == size(torque_history,1));
nominal_max_roll = max(abs(diagnostics.roll_error_deg));
nominal_max_H = max(abs(diagnostics.H));
nominal_max_h_tau = max(abs(diagnostics.h_tau));
assert(nominal_max_roll < 1e-6 && nominal_max_H < 1e-7);
assert(nominal_max_h_tau < 1e-14);

% Compare tightened integration against the production tolerances.
fine_options = odeset('RelTol',1e-12,'AbsTol',1e-14,'MaxStep',2.5e-4);
fine_solution = ode45(odefun, [0 T_max], y0, fine_options);
fine_state = deval(fine_solution,t).';
fine_diagnostics = constraint_diagnostics(t,fine_state,params);
refinement_difference_deg = max(max(abs(rad2deg( ...
    fine_state(:,[1,3,5])-state(:,[1,3,5])))));
assert(refinement_difference_deg < 1e-6);

% Explicit Euler with all six derivatives evaluated at the same old state.
% This isolates integration error without velocity resets or torque deadbands.
steps = [1e-4, 5e-5, 2.5e-5];
euler_roll = zeros(3,1); euler_H = zeros(3,1);
euler_angle_error = zeros(3,1);
for j = 1:numel(steps)
    tt = linspace(0,T_max,round(T_max/steps(j))+1);
    yy = zeros(numel(tt),6); yy(1,:) = y0.';
    for k = 1:numel(tt)-1
        dx = odefun(tt(k), yy(k,:).');
        yy(k+1,:) = yy(k,:) + (tt(k+1)-tt(k))*dx.';
    end
    dd = constraint_diagnostics(tt,yy,params);
    reference = deval(fine_solution,tt).';
    euler_roll(j) = max(abs(dd.roll_error_deg));
    euler_H(j) = max(abs(dd.H));
    euler_angle_error(j) = max(max(abs(rad2deg( ...
        yy(:,[1,3,5])-reference(:,[1,3,5])))));
end
assert(all(diff(euler_roll)<0) && all(diff(euler_angle_error)<0));
convergence = table(steps(:),euler_roll,euler_H,euler_angle_error, ...
    'VariableNames',{'step_s','max_Listing_error_deg','max_H_rad_s', ...
    'max_angle_difference_from_ode45_deg'});
writetable(convergence, fullfile(verify_dir,'integration_convergence.csv'));
disp(convergence);

% A nonzero initial roll offset must survive: there is no Listing projection.
offset_initial = y0; offset_initial(5) = offset_initial(5)+0.02;
offset_solution = ode45(odefun,[0 T_max],offset_initial,fine_options);
offset_state = deval(offset_solution,t).';
offset_diagnostics = constraint_diagnostics(t,offset_state,params);
offset_preservation_error = max(abs(deg2rad( ...
    offset_diagnostics.roll_error_deg)-0.02));
assert(offset_preservation_error < 1e-8);

% Eccentric initial gaze, nonzero Listing-compatible initial angular rates.
eccentric_params = params;
eccentric_params.eta0 = deg2rad([10;-8]);
eccentric_params.target = deg2rad([-12;15]);
ee = [eccentric_params.eta0(1);deg2rad(3); ...
      eccentric_params.eta0(2);deg2rad(-2);0;0];
ee(5) = 2*atan(tan(ee(1)/2)*tan(ee(3)/2));
ee(6) = (sin(ee(3))*ee(2)+sin(ee(1))*ee(4)) ...
    /(1+cos(ee(1))*cos(ee(3)));
eccentric_rhs = @(time,y) full_dynamics(time,y,eccentric_params);
eccentric_solution = ode45(eccentric_rhs,[0 T_max],ee,fine_options);
eccentric_state = deval(eccentric_solution,t).';
eccentric_diagnostics = constraint_diagnostics(t,eccentric_state,eccentric_params);
assert(max(abs(eccentric_diagnostics.roll_error_deg)) < 1e-6);

% Check the matrix-form full ODE against explicit component equations,
% including states that do not satisfy the Listing velocity constraint.
rng(1005,'twister'); component_error = 0;
for j=1:100
    z = randn(6,1); z([1,3,5]) = 0.4*tanh(z([1,3,5]));
    [dx,ttau] = full_dynamics(0.03,z,params);
    p=z(1); th=z(3); dp=z(2); dth=z(4); dph=z(6);
    beta=ttau(1)*cos(p)+ttau(2)*sin(p);
    alpha=-ttau(1)*sin(p)+ttau(2)*cos(p);
    acc=[(beta*tan(th)+ttau(3))/I0 + dth*(dph+sin(th)*dp)/cos(th); ...
         alpha/I0-cos(th)*dp*dph; ...
         beta/(I0*cos(th))+dth*(sin(th)*dph+dp)/cos(th)];
    component_error=max(component_error,max(abs(dx([2,4,6])-acc)));
end
assert(component_error<1e-9);

report.mass_kg=m;
report.radius_m=r;
report.I0_kg_m2=I0;
report.nominal_max_listing_error_deg=nominal_max_roll;
report.nominal_max_H_rad_s=nominal_max_H;
report.nominal_max_h_tau_Nm=nominal_max_h_tau;
report.refined_max_listing_error_deg=max(abs(fine_diagnostics.roll_error_deg));
report.refinement_max_angle_change_deg=refinement_difference_deg;
report.final_angles_deg=rad2deg(state(end,[1,3,5]));
report.final_Listing_roll_deg=rad2deg(diagnostics.listing_roll(end));
report.nonzero_initial_roll_offset_preservation_error_rad=offset_preservation_error;
report.eccentric_max_listing_error_deg=max(abs(eccentric_diagnostics.roll_error_deg));
report.component_equation_max_error_rad_s2=component_error;
report.integration_convergence=table2struct(convergence);
% Reproducibility metadata for the Biomimetics Figure 5 rerun.
% The simulation script and all numerical checks above are unchanged.
report.matlab_version=version;
report.solver='ode45 (Dormand-Prince explicit Runge-Kutta)';
report.relative_tolerance=params.RelTol;
report.absolute_tolerance=params.AbsTol;
report.configured_maximum_internal_step_s=params.MaxStep;
report.maximum_actual_internal_step_s=max(diff(solution.x));
report.minimum_actual_internal_step_s=min(diff(solution.x));
report.accepted_step_count=numel(solution.x)-1;
report.saved_sample_spacing_s=dt;
report.integration_interval_s=[0 T_max];
report.plot_interval_s=[0 comparison_end];
report.saved_sample_count=numel(t);
report.proportional_gains_s_minus_2=[k1_p k2_p];
report.derivative_gains_s_minus_1=[k1_d k2_d];
report.reference_time_constant_s=tau_ref;
report.yaw_pitch_target_deg=rad2deg(params.target).';
report.initial_state=y0.';
report.cartesian_torque_order={'tau_phi (x)','tau_theta (y)','tau_psi (z)'};
report.plot_torque_order={'tau_psi (z)','tau_theta (y)','tau_phi (x)'};
report.peak_absolute_torque_Nm_cartesian=max(abs(torque_history),[],1);
report.simulation_code='code/figure5/Dynamic_model1210.m';
report.verification_code='code/figure5/verify_dynamic_model_integration.m';
report.roll_and_roll_rate_integrated_independently=true;
report.roll_projection_used=false;
report.torque_deadband_used=false;
report.roll_rate_clamp_used=false;
report.all_checks_passed=true;
fid=fopen(fullfile(verify_dir,'integration_checks.json'),'w','n','UTF-8');
fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true)); fclose(fid);
disp(jsonencode(report,'PrettyPrint',true));
end

function restore_figures(previous_visibility, existing_figures)
set(groot, 'defaultFigureVisible', previous_visibility);
new_figures = setdiff(findall(groot, 'Type', 'figure'), existing_figures);
delete(new_figures);
end
