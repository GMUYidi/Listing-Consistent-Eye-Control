function out = fit_retained_segment(channel, fname, model_type)
% Preserved two-stage, 20-start numerical fitter with public CSV input.
% The historical routine computes both branches even when underdamped output
% is requested. That sequence is retained, including its random-number use.
% Multiple starts do not establish a global optimum.
if nargin < 3, model_type = 'underdamp'; end
model_type = lower(model_type);
S = readtable(fname);
t = S.time_s;
y = S.angle_deg;
    t = t(:); y = y(:);
    assert(numel(t)==numel(y) && numel(t)>=5, '时间/信号长度异常。');

    if numel(y)>15
        y_s = movmean(y, max(3, round(0.02*numel(y))));
    else
        y_s = y;
    end
    y_s = y_s(:);

    n = numel(t);
    y0_est   = median(y_s(1:max(3,round(0.15*n))));
    yinf_est = median(y_s(max(1,round(0.85*n)):end));
    dir = sign(yinf_est - y0_est);  % +1: 上升；-1: 下降
    if dir==0, dir = 1; end

    y_fit_sig = dir * y_s;
    y0_fit    = dir * y_s(1);
    yinf_fit  = dir * yinf_est;
    A_fit0    = max(yinf_fit - y0_fit, 1e-6);

    y5  = y0_fit + 0.05*A_fit0;
    i5  = find(y_fit_sig >= y5,  1, 'first'); if isempty(i5), i5 = 1; end
    i98 = find(y_fit_sig >= (y0_fit + 0.98*A_fit0), 1, 'first'); if isempty(i98), i98 = n; end

    i1 = max(1, i5 - 3); i2 = i98;
    if i2 < i1, i1 = 1; i2 = n; end

    t_roi = t(i1:i2);         t_roi = t_roi(:);
    y_roi = y_fit_sig(i1:i2); y_roi = y_roi(:);
    assert(numel(t_roi)==numel(y_roi) && ~isempty(t_roi), 'ROI 长度为空或不一致');

    t0_0 = t(i5);
    Ts0  = max(t(i98) - t0_0, 0.10);

    sslow   = -4/max(Ts0,1e-3);
    sfast   = 3*sslow;
    p0_over = [sslow, sfast, t0_0, y0_fit, A_fit0];
    lb_over = [-Inf, -Inf, t(1)-0.2*range(t),  y0_fit, -Inf];
    ub_over = [ -1e-6,  -1e-6, t(end),         y0_fit,  Inf];

    [p_over, ~] = do_global_fit(@overdamped_model, p0_over, lb_over, ub_over, t_roi, y_roi);

    lb2 = [-Inf, -Inf, t(1)-0.2*range(t)];
    ub2 = [ -1e-6, -1e-6, t(end)];
    yinf_fixed   = mean(y_s(t >= (t(end)-0.10)));   % 原符号域尾段均值
    y0_fixed_fit = p_over(4);                       % = y0_fit
    A_fixed_fit  = dir * (yinf_fixed - (y0_fixed_fit*dir));

    [p_over2, ~] = do_global_fit( ...
        @(tt,pp) overdamped_model(tt, [pp(1) pp(2) pp(3) y0_fixed_fit A_fixed_fit]), ...
        [p_over(1), p_over(2), p_over(3)], lb2, ub2, t_roi, y_roi);

    s_over  = sort(p_over2(1:2));
    s1_over = s_over(1); s2_over = s_over(2);
    t0_over = p_over2(3);
    y0_over = y0_fixed_fit;
    A_over  = A_fixed_fit;

    y_over_fit_full = overdamped_model(t, [s1_over s2_over t0_over y0_over A_over]) * dir;
    rmse_over = sqrt(mean((y_over_fit_full - y_s).^2));
    mae_over  = mean(abs(y_over_fit_full - y_s));
    kd_over = -(s1_over + s2_over);
    kp_over =  (s1_over * s2_over);

    [zeta0, wn0] = guess_underdamped_init(t, y_fit_sig, t0_0, Ts0);
    p0_under   = [zeta0, wn0, t0_0, y0_fit, A_fit0];
    lb_under   = [ 0.02,  1e-3, t(1)-0.2*range(t),  y0_fit, -Inf];
    ub_under   = [ 0.98,  1e4,  t(end),            y0_fit,  Inf];

    [p_under, ~] = do_global_fit(@underdamped_model, p0_under, lb_under, ub_under, t_roi, y_roi);

    [p_under2, ~] = do_global_fit( ...
        @(tt,pp) underdamped_model(tt, [pp(1) pp(2) pp(3) y0_fit A_fit0]), ...
        [p_under(1), p_under(2), p_under(3)], ...
        lb_under(1:3), ub_under(1:3), t_roi, y_roi);

    zeta = p_under2(1); wn = p_under2(2); t0_under = p_under2(3);
    y0_under = y0_fit;  A_under = A_fit0;

    s_real_u = -zeta*wn;
    s_im_u   =  wn*sqrt(max(1 - zeta^2, 0));
    s1_under = s_real_u + 1i*s_im_u;
    s2_under = s_real_u - 1i*s_im_u;

    y_under_fit_full = underdamped_model(t, [zeta wn t0_under y0_under A_under]) * dir;
    rmse_under = sqrt(mean((y_under_fit_full - y_s).^2));
    mae_under  = mean(abs(y_under_fit_full - y_s));
    kd_under = 2*zeta*wn;
    kp_under = wn^2;

    switch model_type
        case 'overdamp'
            chosen   = 'overdamped';
            y_best   = y_over_fit_full;
            kp       = kp_over;
            kd       = kd_over;
        case 'underdamp'
            chosen   = 'underdamped';
            y_best   = y_under_fit_full;
            kp       = kp_under;
            kd       = kd_under;
        otherwise  % 'auto'
            if rmse_over <= rmse_under
                chosen = 'overdamped';
                y_best = y_over_fit_full;
                kp = kp_over; kd = kd_over;
            else
                chosen = 'underdamped';
                y_best = y_under_fit_full;
                kp = kp_under; kd = kd_under;
            end
    end

    out = struct();
    out.channel   = channel;
    out.direction = dir;
    out.model_type = model_type;
    out.over = struct('s1',s1_over,'s2',s2_over,'t0',t0_over,'y0',y0_over*dir,'A',A_over*dir, ...
                      'kp',kp_over,'kd',kd_over,'rmse',rmse_over,'mae',mae_over);
    out.under= struct('zeta',zeta,'wn',wn,'t0',t0_under,'y0',y0_under*dir,'A',A_under*dir, ...
                      'kp',kp_under,'kd',kd_under,'rmse',rmse_under,'mae',mae_under, ...
                      's1',s1_under,'s2',s2_under);
    out.chosen = chosen;   % 最终真正用来取极点的模型
    out.kp = kp; out.kd = kd;
end


function y = overdamped_model(t, p)
    t = t(:);
    s1 = p(1); s2 = p(2); t0 = p(3); y0 = p(4); A = p(5);
    tt = max(t - t0, 0);
    e  = (s2/(s2 - s1)) .* exp(s1*tt) - (s1/(s2 - s1)) .* exp(s2*tt);
    y  = y0 + A * (1 - e);
    y  = y(:);
end

function y = underdamped_model(t, p)
    t = t(:);
    zeta = p(1); wn = p(2); t0 = p(3); y0 = p(4); A = p(5);
    tt = max(t - t0, 0);
    wd  = wn * sqrt(max(1 - zeta^2, 1e-9));
    phi = atan2(sqrt(max(1 - zeta^2, 0)), zeta);
    resp = 1 - (1/sqrt(max(1 - zeta^2,1e-9))) .* exp(-zeta*wn*tt) .* sin(wd*tt + phi);
    y = y0 + A * resp;
    y = y(:);
end

function [p_best, used_lsq_any] = do_global_fit(model_t_p, p0, lb, ub, t, y)

    t = t(:); y = y(:);
    dim = numel(p0);
    Nmulti = 20;  % 可按需加大，比如 30/50

    p_all   = zeros(Nmulti, dim);
    f_all   = inf(Nmulti, 1);
    used_lsq_any = false;

    lb_f = lb; ub_f = ub;
    big = 1e3;
    idx_inf_low  = isinf(lb_f) & lb_f<0;
    idx_inf_high = isinf(ub_f) & ub_f>0;
    lb_f(idx_inf_low)  = -big;
    ub_f(idx_inf_high) =  big;

    for k = 1:Nmulti
        if k == 1
            p_init = p0(:);
        else
            r = rand(dim,1);
            p_init = lb_f(:) + r .* (ub_f(:) - lb_f(:));
        end

        [p_k, used_lsq_k] = do_fit(model_t_p, p_init, lb, ub, t, y);
        used_lsq_any = used_lsq_any | used_lsq_k;
        y_k = model_t_p(t, p_k);
        f_all(k) = sqrt(mean((y_k - y).^2));  % RMSE
        p_all(k,:) = p_k(:).';
    end

    [~, idx_best] = min(f_all);
    p_best = p_all(idx_best,:).';
end

function [p_opt, used_lsq] = do_fit(model_t_p, p0, lb, ub, t, y)
    used_lsq = false;
    t = t(:); y = y(:);
    fun = @(pp, tt) model_t_p(tt, pp);
    use_lsq = license('test','Optimization_Toolbox') || ~isempty(ver('optim'));
    if use_lsq
        opts = optimoptions('lsqcurvefit','Display','none', ...
                            'MaxFunctionEvaluations',5e4,'MaxIterations',2e3);
        p_opt = lsqcurvefit(fun, p0, t, y, lb, ub, opts);
        used_lsq = true;
    else
        obj = @(pp) sqrt(mean((model_t_p(t, proj(pp,lb,ub)) - y).^2));
        p_opt = fminsearch(obj, p0, optimset('Display','off', ...
                          'MaxFunEvals',5e4,'MaxIter',2e3));
        p_opt = proj(p_opt, lb, ub);
    end
end

function [zeta0, wn0] = guess_underdamped_init(t, y_up, t0_0, Ts0)
    t = t(:); y_up = y_up(:);
    Mp = 0;
    try
        [pks, ~] = findpeaks(y_up, t, 'MinPeakProminence', 0.002);
        if ~isempty(pks)
            yinf = median(y_up(max(1,round(0.85*numel(y_up))):end));
            overshoot = max(0, max(pks) - yinf);
            Mp = overshoot / max(abs(yinf - y_up(1)), 1e-9);
        end
    catch
        [maxv, ~] = max(y_up);
        yinf = median(y_up(max(1,round(0.85*numel(y_up))):end));
        overshoot = max(0, maxv - yinf);
        Mp = overshoot / max(abs(yinf - y_up(1)), 1e-9);
    end

    if Mp > 0.005
        zfun = @(z) exp(-z*pi/sqrt(max(1-z^2,1e-9))) - Mp;
        try
            zeta0 = fzero(zfun, [0.05, 0.95]);
        catch
            zeta0 = 0.7;
        end
    else
        zeta0 = 0.7;
    end

    wn0 = 4 / max(zeta0*Ts0, 1e-3);  % Ts ≈ 4/(ζ ωn)
    wn0 = min(max(wn0, 5), 2000);
end

function v = getf(S,name)
    if ~isfield(S,name), error('缺少字段 %s',name); end
    v = S.(name);
end

function s = ternary(cond,a,b)
    if cond, s=a; else, s=b; end
end

function p = proj(p, lb, ub)
    p = max(p, lb); p = min(p, ub);
end
