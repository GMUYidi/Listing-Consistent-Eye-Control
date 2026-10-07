function render_figure5
% Same two-panel plotting block as the corrected source, using saved states.
% Only change: torque y-limit [-0.3,0.2] to avoid clipping the negative peak.
output_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'outputs', 'figure5');
load(fullfile(output_dir,'feedback_result_corrected.mat'), ...
    't','x1','x3','x5','tau_x_values','tau_y_values','tau_z_values','params');
x1_deg=rad2deg(x1); x3_deg=rad2deg(x3); x5_deg=rad2deg(x5);
x1_d_deg=rad2deg(params.target(1)); x3_d_deg=rad2deg(params.target(2));
x5_d_deg=rad2deg(2*atan(tan(params.target(1)/2)*tan(params.target(2)/2)));
previous_visibility=get(groot,'defaultFigureVisible');
cleanup=onCleanup(@() set(groot,'defaultFigureVisible',previous_visibility)); %#ok<NASGU>
set(groot,'defaultFigureVisible','off');
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
% A small y-margin keeps the negative peak inside the axes.
ylim(ax_torque, [-0.3, 0.2]);

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

close(comparison_fig);
end
