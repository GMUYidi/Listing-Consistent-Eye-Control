% Torque-driven, independently integrated full six-state eye dynamics.
% State: [psi; psi_dot; theta; theta_dot; phi; phi_dot].
% Spatial torque order: [tau_phi; tau_theta; tau_psi] = [tau_x; tau_y; tau_z].
% Listing roll/rate are diagnostics only: never projected onto the state.
% Run verify_dynamic_model_integration.m for numerical convergence checks.
clear;
clc;
% Initial state and parameters

x1_initial = deg2rad(0);
x2_initial = 0.00;
x3_initial = deg2rad(0);
x4_initial = 0.00;
x5_initial = asin((sin(x1_initial) * sin(x3_initial)) / (1 + cos(x1_initial) * cos(x3_initial)));
% Consistent initial rate; remains zero for the default rest condition.
x6_initial = (sin(x3_initial)*x2_initial + sin(x1_initial)*x4_initial) ...
    / (1 + cos(x1_initial)*cos(x3_initial));
y0 = [x1_initial; x2_initial; x3_initial; x4_initial; x5_initial; x6_initial];

OLR = [ .0055; -.018; 0];
OMR = [-.0055; -.018; 0];
OSR = [ 0; -.018; .007];
OIR = [ 0; -.018; -.007];
OSO = [-.024; .0125; .007];
OIO = [-.024; .0125; -.007];
% Anatomical coordinates above use +y as gaze; dynamics below use +x.
R_to_gaze_x = [0, 1, 0; -1, 0, 0; 0, 0, 1];
Pull = R_to_gaze_x * [OLR, OMR, OSR, OIR, OSO, OIO];

% Insertion points (3×6)
ILR = [ .0134; .0048;  0];
IMR = [-.0134; .0048; 0];
ISR = [ 0; .0052; .0133];
IIR = [ 0; .0016; -.01425];
ISO = [ 0;  -.0052;   .0133];
IIO = [ 0;   -.0052; -.0133];
Ins_init = R_to_gaze_x * [ILR, IMR, ISR, IIR, ISO, IIO];
EyeBall_Radius = 0.015;

m = 0.0075;        % Platform design mass: 7.5 g (kg)
r = EyeBall_Radius;        % Radius (m)
I0 = 2*m*r^2/5;    % Homogeneous spherical inertia (kg m^2)
% [k1_p, k1_d] = pd_from_spec(0.02, 0.05); %pd_from_spec(T_r, T_s)
% [k2_p, k2_d] = pd_from_spec(0.02, 0.05); %pd_from_spec(T_r, T_s)
% s1_yaw = -62.7;
% s2_yaw = -62.7;
% s1_pitch = -47;
% s2_pitch = -46.9;
% [k1_p, k1_d] = pd_from_poles(s1_yaw, s2_yaw); %pd_from_poles(s1, s2)
% [k2_p, k2_d] = pd_from_poles(s1_pitch, s2_pitch); %pd_from_poles(s1, s2)

s1 = -50+35i;
s2 = -50-35i;
s3 = -50+40i;
s4 = -50-40i;
[k1_p, k1_d] = pd_from_poles(s1, s2); %pd_from_poles(s1, s2)
[k2_p, k2_d] = pd_from_poles(s3, s4); %pd_from_poles(s1, s2)
% Target state
x1_d = deg2rad(-28.7);     % Final yaw target
x3_d = deg2rad(10.8);      % Final pitch target
x5_d = asin((sin(x1_d) * sin(x3_d)) / (1 + cos(x1_d) * cos(x3_d))); % Target angle for x5

tau_ref = 0.02;  % 20 ms smoothing time constant


% dt controls saved/plot sample spacing, not ode45's internal integration step.
dt = 0.0001;
T_max = 0.5;
output_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'outputs', 'figure5');
if ~isfolder(output_dir), mkdir(output_dir); end
params.I0 = I0;
params.KP = diag([k1_p, k2_p]);
params.KD = diag([k1_d, k2_d]);
params.eta0 = y0([1,3]);
params.target = [x1_d; x3_d];
params.tau_ref = tau_ref;
params.RelTol = 1e-10;
params.AbsTol = 1e-12;
params.MaxStep = 1e-3;
params.T_max = T_max;

% Handles are retained in the workspace for the companion verification script.
full_dynamics = @full_eye_rhs;
constraint_diagnostics = @eye_constraint_diagnostics;
odefun = @(time, state) full_dynamics(time, state, params);
ode_options = odeset('RelTol', params.RelTol, 'AbsTol', params.AbsTol, ...
    'MaxStep', params.MaxStep);
time_samples = unique([0:dt:T_max, T_max]);
solution = ode45(odefun, [0 T_max], y0, ode_options);
t = time_samples;
state = deval(solution, t).';
x1 = state(:,1).'; x2 = state(:,2).';
x3 = state(:,3).'; x4 = state(:,4).';
x5 = state(:,5).'; x6 = state(:,6).';

% All histories use the same sample times, including the final sample.
diagnostics = constraint_diagnostics(t, state, params);
torque_history = diagnostics.torque;
tau_x_values = torque_history(:,1).';
tau_y_values = torque_history(:,2).';
tau_z_values = torque_history(:,3).';
x2dot_values = diagnostics.acceleration(:,1).';
x4dot_values = diagnostics.acceleration(:,2).';
x6dot_values = diagnostics.acceleration(:,3).';
x6_gt_values = diagnostics.listing_rate.';
x6_err_values = x6 - x6_gt_values;

fprintf('\nFull-state integration, m = %.4g kg, r = %.4g m\n', m, r);
fprintf('Max |phi - L(psi,theta)|: %.6e deg\n', ...
    max(abs(diagnostics.roll_error_deg)));
fprintf('Max |H|:                  %.6e rad/s\n', max(abs(diagnostics.H)));
fprintf('Max |h^T tau|:            %.6e N m\n', max(abs(diagnostics.h_tau)));

% Preserve the legacy force allocator as an optional, separate diagnostic.
% It is not used as the input to this torque-level consistency experiment.
% Its geometry/force approximation requires its own actuation validation.
run_force_diagnostic = false;
Force_history = nan(6, numel(t));
allocation_error = nan(3, numel(t));
if run_force_diagnostic
    for k = 1:numel(t)
        Ins_rotated = RotateInsertionPoints(x1(k), x3(k), x5(k), Ins_init);
        [F, Torque_Dir, ~] = Force_Distributor( ...
            torque_history(k,:).', Pull, Ins_rotated, EyeBall_Radius);
        Force_history(:,k) = F(:);
        allocation_error(:,k) = Torque_Dir*(F(:)*EyeBall_Radius) ...
            - torque_history(k,:).';
    end
end

save(fullfile(output_dir, 'feedback_result_corrected.mat'), 'params', ...
    'm', 'r', 'I0', 't', 'state', 'torque_history', 'diagnostics', ...
    'x1', 'x2', 'x3', 'x4', 'x5', 'x6', ...
    'tau_x_values', 'tau_y_values', 'tau_z_values', ...
    'k1_p', 'k1_d', 'k2_p', 'k2_d', 'tau_ref');
state_table = array2table([t(:), state], 'VariableNames', ...
    {'time_s','psi_rad','psi_dot_rad_s','theta_rad','theta_dot_rad_s', ...
    'phi_rad','phi_dot_rad_s'});
writetable(state_table, fullfile(output_dir, 'feedback_states_corrected.csv'));
torque_table = array2table([t(:), torque_history], 'VariableNames', ...
    {'time_s','tau_phi_x_Nm','tau_theta_y_Nm','tau_psi_z_Nm'});
writetable(torque_table, fullfile(output_dir, 'feedback_torques_corrected.csv'));
residual_table = table(t(:), diagnostics.roll_error_deg, diagnostics.H, ...
    diagnostics.h_tau, 'VariableNames', ...
    {'time_s','Listing_roll_error_deg','H_rad_s','h_tau_Nm'});
writetable(residual_table, fullfile(output_dir, 'feedback_residuals_corrected.csv'));

x5_L = asin((sin(x1) .* sin(x3)) ./ (1 + cos(x1) .* cos(x3)));
% Plot the evolution of x1, x3, x5 over time
x1_deg = rad2deg(x1);
x3_deg = rad2deg(x3);
x5_deg = rad2deg(x5);
x5_L_deg = rad2deg(x5_L);
x1_d_deg = rad2deg(x1_d);
x3_d_deg = rad2deg(x3_d);
x5_d_deg = rad2deg(x5_d);

% Plot angles in degrees
fontsize = 18;
linewidth = 3;

figure;
plot(t, x1_deg, 'r-', 'LineWidth', linewidth); hold on;
plot(t, x3_deg, 'b-', 'LineWidth', linewidth);
plot(t, x5_deg, 'g-', 'LineWidth', linewidth);
plot(t, x5_L_deg, 'y-', 'LineWidth', linewidth);
h1 = plot(t, ones(size(t)) * x1_d_deg, 'r:', 'LineWidth', linewidth);
h2 = plot(t, ones(size(t)) * x3_d_deg, 'b:', 'LineWidth', linewidth);
h3 = plot(t, ones(size(t)) * x5_d_deg, 'g:', 'LineWidth', linewidth);

% ===== 在目标线末端标数字 =====
% 选择末端点的位置（比如 t(end)）
xpos = t(end) + 0.01;   % 稍微往右偏一点避免重叠

text(xpos, x1_d_deg, sprintf('%.1f°', x1_d_deg), ...
    'Color', 'r', 'FontSize', fontsize, 'HorizontalAlignment','left');

text(xpos, x3_d_deg, sprintf('%.1f°', x3_d_deg), ...
    'Color', 'b', 'FontSize', fontsize, 'HorizontalAlignment','left');

text(xpos, x5_d_deg, sprintf('%.1f°', x5_d_deg), ...
    'Color', 'g', 'FontSize', fontsize, 'HorizontalAlignment','left');

% Figure 5 is a simulation-only illustration; experimental data are separate.

xlabel('Time (s)', 'FontSize', fontsize);
ylabel('Angle (deg)', 'FontSize', fontsize);
grid on;
set(gca, 'FontSize', fontsize);


% Diagnostic generalized forces: Q = J_omega(psi,theta)' * tau_world.
psi_plot = x1;
theta_plot = x3;
Q_psi_values = tau_z_values;
Q_theta_values = -sin(psi_plot).*tau_x_values + cos(psi_plot).*tau_y_values;
Q_phi_values = cos(psi_plot).*cos(theta_plot).*tau_x_values + ...
    sin(psi_plot).*cos(theta_plot).*tau_y_values - sin(theta_plot).*tau_z_values;

torque_fig = figure('Position', [100, 100, 800, 420]);
h_Q_psi = plot(t, Q_psi_values, 'r-', 'LineWidth', linewidth); hold on;
h_Q_theta = plot(t, Q_theta_values, 'b-', 'LineWidth', linewidth);
h_Q_phi = plot(t, Q_phi_values, 'g-', 'LineWidth', linewidth);
legend([h_Q_psi, h_Q_theta, h_Q_phi], ...
    'Q_\psi', 'Q_\theta', 'Q_\phi', ...
    'FontSize', 14, 'Location', 'northeast');
xlabel('Time (s)', 'FontSize', 16);
ylabel('Generalized torque (N·m)', 'FontSize', 16);
xlim([0, min(0.3, t(end))]);
grid on;
set(gca, 'FontSize', 14);

% In the paper's fixed-axis convention, tau_psi, tau_theta, tau_phi are
% the z-, y-, and x-axis inputs to the spherical dynamics, respectively.
comparison_fig = figure('Units', 'inches', 'Position', [1, 1, 7.5, 5.5], 'Color', 'w');
comparison_layout = tiledlayout(comparison_fig, 2, 1, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
paper_curve_width = 3.0;
paper_target_width = 1.8;
color_yaw = [0.82, 0.18, 0.20];
color_pitch = [0.13, 0.32, 0.76];
color_torsion = [0.12, 0.58, 0.32];
comparison_end = min(0.3, t(end));

ax_torque = nexttile(comparison_layout);
plot(ax_torque, t, 1e3*tau_z_values, 'Color', color_yaw, ...
    'LineWidth', paper_curve_width); hold(ax_torque, 'on');
plot(ax_torque, t, 1e3*tau_y_values, 'Color', color_pitch, ...
    'LineWidth', paper_curve_width);
plot(ax_torque, t, 1e3*tau_x_values, 'Color', color_torsion, ...
    'LineWidth', paper_curve_width);
title(ax_torque, '(a) Net torque (fixed-axis convention)', ...
    'FontWeight', 'normal', 'FontSize', 18);
ylabel(ax_torque, 'Torque (mN\cdotm)', 'FontSize', 17);
legend(ax_torque, {'\tau_\psi', '\tau_\theta', '\tau_\phi'}, ...
    'Location', 'northeast', 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 15);
grid(ax_torque, 'on');

ax_orientation = nexttile(comparison_layout);
target_yaw_color = 0.55*[1, 1, 1] + 0.45*color_yaw;
target_pitch_color = 0.55*[1, 1, 1] + 0.45*color_pitch;
target_torsion_color = 0.55*[1, 1, 1] + 0.45*color_torsion;
plot(ax_orientation, [0, comparison_end], [x1_d_deg, x1_d_deg], '--', ...
    'Color', target_yaw_color, 'LineWidth', paper_target_width, 'HandleVisibility', 'off');
hold(ax_orientation, 'on');
plot(ax_orientation, [0, comparison_end], [x3_d_deg, x3_d_deg], '--', ...
    'Color', target_pitch_color, 'LineWidth', paper_target_width, 'HandleVisibility', 'off');
plot(ax_orientation, [0, comparison_end], [x5_d_deg, x5_d_deg], '--', ...
    'Color', target_torsion_color, 'LineWidth', paper_target_width, 'HandleVisibility', 'off');
h_psi_orientation = plot(ax_orientation, t, x1_deg, 'Color', color_yaw, ...
    'LineWidth', paper_curve_width);
h_theta_orientation = plot(ax_orientation, t, x3_deg, 'Color', color_pitch, ...
    'LineWidth', paper_curve_width);
h_phi_orientation = plot(ax_orientation, t, x5_deg, 'Color', color_torsion, ...
    'LineWidth', paper_curve_width);
h_target_key = plot(ax_orientation, NaN, NaN, '--', ...
    'Color', [0.5, 0.5, 0.5], 'LineWidth', paper_target_width);
title(ax_orientation, '(b) Simulated eye orientation', ...
    'FontWeight', 'normal', 'FontSize', 18);
ylabel(ax_orientation, 'Angle (deg)', 'FontSize', 17);
xlabel(ax_orientation, 'Time (s)', 'FontSize', 17);
legend(ax_orientation, [h_psi_orientation, h_theta_orientation, h_phi_orientation, h_target_key], ...
    {'\psi (yaw)', '\theta (pitch)', '\phi (torsion)', 'Target (dashed)'}, ...
    'Location', 'northwest', 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 15);
grid(ax_orientation, 'on');

linkaxes([ax_torque, ax_orientation], 'x');
xlim(ax_torque, [0, comparison_end]);
angle_limits = [x1_deg, x3_deg, x5_deg, x1_d_deg, x3_d_deg, x5_d_deg];
angle_min = min(angle_limits);
angle_max = max(angle_limits);
angle_span = max(angle_max - angle_min, 1);
ylim(ax_orientation, [angle_min - 0.10*angle_span, angle_max + 0.35*angle_span]);
set(ax_torque, 'XTickLabel', []);
set([ax_torque, ax_orientation], 'FontSize', 16, 'LineWidth', 1.0, 'FontName', 'Arial');
exportgraphics(comparison_fig, fullfile(output_dir, 'net_torque_eye_orientation.png'), 'Resolution', 600);
exportgraphics(comparison_fig, fullfile(output_dir, 'net_torque_eye_orientation.pdf'), 'ContentType', 'vector');

% figure;
% plot(t, x2dot_values, 'r-', 'LineWidth', linewidth); hold on;
% plot(t, x4dot_values, 'b-', 'LineWidth', linewidth);
% plot(t, x6dot_values, 'm-', 'LineWidth', linewidth);
% legend('\dot{x}_2', '\dot{x}_4', '\dot{x}_6', 'Location', 'best', 'FontSize', fontsize);
% xlabel('Time (s)', 'FontSize', fontsize);
% ylabel('Angular acceleration (rad/s²)', 'FontSize', fontsize);
% title('Time evolution of \dot{x}_2, \dot{x}_4, \dot{x}_6', 'FontSize', fontsize);
% grid on;
% set(gca, 'FontSize', fontsize);
%
% % Convert angular velocity from rad/s to deg/s
% x2_deg = rad2deg(x2);
% x4_deg = rad2deg(x4);
% x6_deg = rad2deg(x6);

% figure;
% plot(t, x2_deg, 'r-', 'LineWidth', linewidth); hold on;
% plot(t, x4_deg, 'b-', 'LineWidth', linewidth);
% plot(t, x6_deg, 'm-', 'LineWidth', linewidth);
% legend('x_2 (deg/s)', 'x_4 (deg/s)', 'x_6 (deg/s)', 'Location', 'best', 'FontSize', fontsize);
% xlabel('Time (s)', 'FontSize', fontsize);
% ylabel('Angular velocity (deg/s)', 'FontSize', fontsize);
% title('Time evolution of x_2, x_4, x_6', 'FontSize', fontsize);
% grid on;
% set(gca, 'FontSize', fontsize);

if run_force_diagnostic
figure;
plot(t, Force_history', 'LineWidth', linewidth);
legend('F_{LR}','F_{MR}','F_{SR}','F_{IR}','F_{SO}','F_{IO}', 'Location','best', 'FontSize', fontsize);
xlabel('Time (s)', 'FontSize', fontsize);
ylabel('Force (N)', 'FontSize', fontsize);
title('Time evolution of muscle forces', 'FontSize', fontsize);
grid on;
set(gca, 'FontSize', fontsize);
end
%%%%%%%%%%%%%%%%
% figure;
%
% % 计算三组差值
% diff_LR_MR = Force_history(1,:) - Force_history(2,:);
% diff_SR_IR = Force_history(3,:) - Force_history(4,:);
% diff_SO_IO = Force_history(5,:) - Force_history(6,:);
%
% % 绘图
% plot(t, diff_LR_MR, 'LineWidth', linewidth); hold on;
% plot(t, diff_SR_IR, 'LineWidth', linewidth);
% plot(t, diff_SO_IO, 'LineWidth', linewidth);
%
% legend('F_{LR} - F_{MR}', 'F_{SR} - F_{IR}', 'F_{SO} - F_{IO}', ...
%        'Location','best', 'FontSize', fontsize);
%
% xlabel('Time (s)', 'FontSize', fontsize);
% ylabel('Force Difference (N)', 'FontSize', fontsize);
% title('Time Evolution of Muscle Force Differences', 'FontSize', fontsize);
% grid on;
% set(gca, 'FontSize', fontsize);
%%%%%%%%%%%%%%%%
% ==== 构造平滑参考轨迹 ====
alpha = 1 - exp(-t / tau_ref);   % 和仿真里用的一样的 smooth 系数

x1_d_traj = y0(1) + (x1_d-y0(1)) * alpha;        % x1_d_now 随时间的轨迹
x3_d_traj = y0(3) + (x3_d-y0(3)) * alpha;        % x3_d_now 随时间的轨迹

% 对应的 Listing-law 目标 x5_d_now 轨迹
x5_d_now_traj = asin( (sin(x1_d_traj) .* sin(x3_d_traj)) ./ ...
                      (1 + cos(x1_d_traj) .* cos(x3_d_traj)) );
figure; hold on; box on;

% 动态目标
plot(t, x1_d_traj,      'r-', 'LineWidth', 2);  % x1_d_now
plot(t, x3_d_traj,      'b-', 'LineWidth', 2);  % x3_d_now
plot(t, x5_d_now_traj,  'g-', 'LineWidth', 2);  % x5_d_now

% 恒定目标
plot(t, x1_d * ones(size(t)), 'r--', 'LineWidth', 2);
plot(t, x3_d * ones(size(t)), 'b--', 'LineWidth', 2);
plot(t, x5_d * ones(size(t)), 'g--', 'LineWidth', 2);

xlabel('Time (s)');
ylabel('Angle (rad)');
title('Dynamic vs Constant Reference Angles');
legend('x_1^d(t)', 'x_3^d(t)', 'x_5^d(t)', ...
       'x_1^{target}',    'x_3^{target}',    'x_5^{target}', ...
       'Location', 'southeast');
grid on;
set(gca, 'FontSize', 16);

%%%%%%%%%%%%%%%%%%%%%




function [dx, tau, omega] = full_eye_rhs(time, state, p)
    % Independent six-state spherical dynamics: I0*omega_dot = tau.
    % Neither phi nor phi_dot is reconstructed or reset here.
    psi = state(1); theta = state(3);
    rates = state([2,4,6]);
    if abs(psi) >= pi/2 || abs(theta) >= pi/2
        error('DynamicEye:Domain', 'Trajectory left the stated Euler/Listing domain.');
    end
    tau = eye_torque_command(time, state, p);
    E = [0, -sin(psi), cos(psi)*cos(theta); ...
         0,  cos(psi), sin(psi)*cos(theta); ...
         1,         0,          -sin(theta)];
    E_psi = [0, -cos(psi), -sin(psi)*cos(theta); ...
             0, -sin(psi),  cos(psi)*cos(theta); ...
             0,         0,                     0];
    E_theta = [0, 0, -cos(psi)*sin(theta); ...
               0, 0, -sin(psi)*sin(theta); ...
               0, 0,                -cos(theta)];
    E_dot = E_psi*rates(1) + E_theta*rates(2);
    acceleration = E \ (tau/p.I0 - E_dot*rates);
    dx = [rates(1); acceleration(1); rates(2); acceleration(2); ...
          rates(3); acceleration(3)];
    omega = E*rates;
end

function tau = eye_torque_command(time, state, p)
    % Original filtered-target PD command, without reference derivatives.
    psi = state(1); psi_dot = state(2);
    theta = state(3); theta_dot = state(4);
    eta = state([1,3]); eta_dot = state([2,4]);
    eta_d = p.eta0 + (p.target-p.eta0)*(1-exp(-time/p.tau_ref));
    v = -p.KP*(eta-eta_d) - p.KD*eta_dot;
    % Explicit reduced acceleration terms, as in the dissertation.
    a = [(theta_dot^2*sin(psi) ...
          +theta_dot*psi_dot*sin(theta)*(2+cos(psi)*cos(theta))) ...
          /((1+cos(psi)*cos(theta))*cos(theta)); ...
         -psi_dot*cos(theta)*(theta_dot*sin(psi)+psi_dot*sin(theta)) ...
          /(1+cos(psi)*cos(theta))];
    B = [(cos(psi)+cos(theta))/(p.I0*(1+cos(psi)*cos(theta))*cos(theta)), ...
         sin(psi)*sin(theta)/(p.I0*(1+cos(psi)*cos(theta))*cos(theta)); ...
         -sin(psi)*sin(theta)/(p.I0*(1+cos(psi)*cos(theta))), ...
         (cos(psi)+cos(theta))/(p.I0*(1+cos(psi)*cos(theta)))];
    u = B \ (v-a); % [tau_psi; tau_theta]
    tau_phi = (sin(theta)*u(1)-cos(theta)*sin(psi)*u(2)) ...
        /(1+cos(psi)*cos(theta));
    tau = [tau_phi; u(2); u(1)];
    % No independent component deadbands: they would break h'*tau = 0.
end

function d = eye_constraint_diagnostics(times, states, p)
    n = numel(times);
    d.torque = zeros(n,3); d.acceleration = zeros(n,3);
    d.H = zeros(n,1); d.h_tau = zeros(n,1);
    psi = states(:,1); theta = states(:,3);
    d.listing_roll = 2*atan(tan(psi/2).*tan(theta/2));
    d.roll_error_deg = rad2deg(states(:,5)-d.listing_roll);
    d.listing_rate = (sin(theta).*states(:,2)+sin(psi).*states(:,4)) ...
        ./ (1+cos(psi).*cos(theta));
    for k = 1:n
        [dx, tau, omega] = full_eye_rhs(times(k), states(k,:).', p);
        h = [1+cos(psi(k))*cos(theta(k)); ...
             sin(psi(k))*cos(theta(k)); -sin(theta(k))];
        d.torque(k,:) = tau.';
        d.acceleration(k,:) = dx([2,4,6]).';
        d.H(k) = h.'*omega;
        d.h_tau(k) = h.'*tau;
    end
end

function [Force_Vector, Torque_Dir, Force_Dir] = Force_Distributor(Desired_Torque, Pull, Ins, EyeBall_Radius)
    Yes_No = zeros(1, size(Pull, 2));
    origin = zeros(3,1);

    for i = 1:size(Pull, 2)
        A = Pull(:, i);
        B = Ins(:, i);
        AB = B - A;
        AC = origin - A;

        t = dot(AC, AB) / dot(AB, AB);
        Yes_No(i) = (t >= 0 && t <= 1);
    end

    %disp("Yes_No is: ");
    %disp(Yes_No);

    % Direction vectors from insertion to pulley
    Ins_to_Pulley = Pull - Ins;
    Ins_to_Pulley_Dir = Ins_to_Pulley ./ vecnorm(Ins_to_Pulley);

    % Torque direction vectors
    Torque_before = cross(Ins, Ins_to_Pulley);
    Torque_Dir = Torque_before ./ vecnorm(Torque_before);

    Norm_Torque_Dir = vecnorm(Torque_Dir);

    % Force direction vectors
    Force_Dir_1 = cross(Torque_before, Ins);
    Force_Dir = Force_Dir_1 ./ vecnorm(Force_Dir_1);

    Norm_Force_Dir = vecnorm(Force_Dir);

    for i = 1:length(Yes_No)
        if ~Yes_No(i)
            Force_Dir(:, i) = Ins_to_Pulley_Dir(:, i);
        end
    end

    %disp("The norm of the Force Direction ");
    %disp(Norm_Force_Dir);

    Ins_Dir = Ins ./ vecnorm(Ins);
    %disp("The norm of the Ins Direction ");
    %disp(vecnorm(Ins_Dir));

    cosines = zeros(1, size(Force_Dir, 2));
    for i = 1:size(Force_Dir, 2)
        a = Force_Dir(:, i);
        b = Ins_Dir(:, i);
        dot_val = dot(a, b);
        norm_a = norm(a);
        norm_b = norm(b);
        if norm_a ~= 0 && norm_b ~= 0
            cosines(i) = dot_val / (norm_a * norm_b);
        else
            cosines(i) = 0;
        end
    end

    %disp("Cosines vector is:");
    %disp(cosines);

    %disp("This is the direction of the Torque");
    %disp(Torque_Dir);

    %disp("This is the direction of the Norm_of_Torque");
    %disp(Norm_Torque_Dir);

    Force_Vector = zeros(1, 6);

    for i = 1:500
        Image_of_Torque = Torque_Dir' * Desired_Torque;
        Norm_Desired_Torque = norm(Desired_Torque);
        max_val = max(Image_of_Torque);
        %indices = find(Image_of_Torque == max_val);
        [~, index] = max(Image_of_Torque);

        % update Force_Vector
        %Force_Vector(indices) = Force_Vector(indices) + max_val ./ (EyeBall_Radius * sqrt(1 - cosines(indices).^2));

        %disp("The force that should be updated:");
        %disp(max_val / EyeBall_Radius);

        % update Desired_Torque
        %Desired_Torque = Desired_Torque - Torque_Dir(:, indices) * max_val;
        Force_Vector(index) = Force_Vector(index) + max_val / (EyeBall_Radius * sqrt(1 - cosines(index)^2));
        Desired_Torque = Desired_Torque - Torque_Dir(:, index) * max_val;
    end
end

function Ins_rotated = RotateInsertionPoints(psi, theta, phi, Ins_init)
    Rz = [cos(psi), -sin(psi), 0;
          sin(psi),  cos(psi), 0;
                0,        0,   1];

    Ry = [ cos(theta), 0, sin(theta);
                   0, 1,          0;
          -sin(theta), 0, cos(theta)];

    Rx = [1,         0,          0;
          0, cos(phi), -sin(phi);
          0, sin(phi),  cos(phi)];

    R = Rz * Ry * Rx;
    Ins_rotated = R * Ins_init;
end

function [k_p, k_d] = pd_from_spec(T_r, T_s)
    % Compute natural frequency and damping ratio
    omega_n = 1.8 / T_r;
    zeta = -log(0.02) / (omega_n * T_s);

    % Compute PD gains
    k_p = omega_n^2;
    k_d = 2 * zeta * omega_n;

    % Print computed values
    fprintf('--- PD Controller Design ---\n');
    fprintf('Desired Rise Time (T_r): %.4f s\n', T_r);
    fprintf('Desired Settling Time (T_s): %.4f s\n', T_s);
    fprintf('Natural Frequency (ω_n): %.4f rad/s\n', omega_n);
    fprintf('Damping Ratio (ζ): %.4f\n', zeta);

    % Classify damping condition
    if zeta < 1
        fprintf('System type: UNDERDAMPED (oscillatory response)\n');
    elseif abs(zeta - 1) < 1e-3
        fprintf('System type: CRITICALLY DAMPED (fastest non-oscillatory)\n');
    else
        fprintf('System type: OVERDAMPED (slow, non-oscillatory)\n');
    end

    fprintf('Computed k_p: %.4f\n', k_p);
    fprintf('Computed k_d: %.4f\n', k_d);
    fprintf('-----------------------------\n\n');
end


function [k_p, k_d] = pd_from_zeta(T_r, zeta)
    omega_n = 1.8 / T_r;
    k_p = omega_n^2;
    k_d = 2 * zeta * omega_n;

    fprintf('Custom damping design: ζ = %.2f\n', zeta);
    fprintf('ω_n = %.2f → k_p = %.2f, k_d = %.2f\n', omega_n, k_p, k_d);
end


function [k_p, k_d] = pd_from_poles(s1, s2, J, B, K0)
%PD_FROM_POLES  Calculate PD controller gains from desired poles
%   For a 2nd-order system: J*y'' + B*y' + K0*y = u,
%   with PD control: u = k_p*e + k_d*e'
%
%   Inputs:
%       s1, s2 : desired closed-loop poles (can be real or complex)
%       J  : inertia (default 1)
%       B  : viscous damping (default 0)
%       K0 : passive stiffness (default 0)
%
%   Outputs:
%       k_p : proportional gain
%       k_d : derivative gain
%
%   Note:
%       - For real distinct poles (overdamped): s1,s2<0
%       - For complex conjugate poles (underdamped): s2 = conj(s1)

    % Default values
    if nargin < 3 || isempty(J),  J  = 1;  end
    if nargin < 4 || isempty(B),  B  = 0;  end
    if nargin < 5 || isempty(K0), K0 = 0;  end

    % Check stability
    if real(s1) >= 0 || real(s2) >= 0
        error('Poles must be in the left-half plane (negative real part).');
    end

    % Case 1: complex conjugate poles
    if ~isreal(s1) || ~isreal(s2)
        if abs(s2 - conj(s1)) > 1e-9
            error('Complex poles must be conjugate.');
        end
        omega_n = abs(s1);
        zeta = -real(s1)/omega_n;

    % Case 2: real poles (distinct or equal)
    else
        p1 = -s1;  % positive decay rate
        p2 = -s2;
        omega_n = sqrt(p1*p2);
        zeta = (p1 + p2) / (2*omega_n);
    end

    % Compute PD gains
    k_d = 2*J*zeta*omega_n - B;
    k_p =     J*omega_n^2  - K0;
end
