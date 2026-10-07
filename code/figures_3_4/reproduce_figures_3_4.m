function reproduce_figures_3_4(output_dir, replot_only)
% Replot the selected trajectories and perform a deterministic repository refit.
% This does not claim to recover the unrecorded historical random state.
code_dir = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(code_dir));
if nargin < 1, output_dir = fullfile(root,'outputs','figures_3_4'); end
if nargin < 2, replot_only = false; end
if ~isfolder(output_dir), mkdir(output_dir); end
addpath(code_dir);
old_visibility = get(groot,'defaultFigureVisible');
restore_visibility = onCleanup(@() set(groot,'defaultFigureVisible',old_visibility)); %#ok<NASGU>
set(groot,'defaultFigureVisible','off');
segments = readtable(fullfile(root,'data','fitting_segments.csv'),'TextType','string');
assert(height(segments)==17 && sum(segments.channel=="yaw")==9 && sum(segments.channel=="pitch")==8);
if replot_only
    fit_table=readtable(fullfile(output_dir,'repository_refit_parameters.csv'),'TextType','string');
else
seed = 20261007; % Repository rerun seed, not the historical fitting seed.
rng(seed,'twister');
config = struct('purpose','Deterministic repository rerun; not historical random-state recovery', ...
    'seed',seed,'generator','twister','matlab_version',version, ...
    'matlab_release',version('-release'), ...
    'optimization_toolbox_available',license('test','Optimization_Toolbox') || ~isempty(ver('optim')), ...
    'findpeaks_available',exist('findpeaks','file')~=0, ...
    'selected_branch','underdamp','random_starts_per_stage',20);
write_json(fullfile(output_dir,'run_configuration.json'),config);

rows = struct([]);
for k=1:height(segments)
    sid=segments.segment_id(k); axis_name=segments.channel(k);
    data_file=fullfile(root,segments.segment_csv(k));
    fprintf('Refitting %d/%d: %s\n',k,height(segments),sid);
    out=fit_retained_segment(char(axis_name),char(data_file),'underdamp');
    p=out.under;
    rows(k).segment_id=char(sid);
    rows(k).participant_id=char(segments.participant_id(k));
    rows(k).channel=char(axis_name);
    rows(k).zeta=p.zeta;
    rows(k).omega_n_per_s=p.wn;
    rows(k).onset_shift_s=p.t0;
    rows(k).baseline_deg=p.y0;
    rows(k).amplitude_deg=p.A;
    rows(k).pole_real_per_s=real(p.s1);
    rows(k).pole_imag_abs_per_s=abs(imag(p.s1));
    rows(k).kp_per_s2=p.kp;
    rows(k).kd_per_s=p.kd;
    rows(k).rmse_all_filtered_samples_deg=p.rmse;
    rows(k).mae_all_filtered_samples_deg=p.mae;
end
fit_table=struct2table(rows);
writetable(fit_table,fullfile(output_dir,'repository_refit_parameters.csv'));
end
for axis_name=["yaw","pitch"]
    plot_selected_trajectories(root,segments(segments.channel==axis_name,:),axis_name,output_dir);
    plot_response_poles(fit_table(string(fit_table.channel)==axis_name,:),axis_name,output_dir);
end
fprintf('Completed repository refit of %d segments. Outputs: %s\n',height(segments),output_dir);
end

function plot_selected_trajectories(root,meta,axis_name,output_dir)
% Display processing reconstructed from the retained plotting scripts.
% It is deliberately separate from the three-sample fitting preprocessing.
if axis_name=="yaw", end_ms=180; initial_limit=5;
else, end_ms=140; initial_limit=10; end
dense_ms=(0:end_ms)'; dense_s=dense_ms/1000;
fig=figure('Color','w','Position',[100 100 900 600]);
ax=axes(fig); hold(ax,'on'); grid(ax,'on'); colors=hsv(height(meta));
display=table(dense_ms,'VariableNames',{'time_ms'});
for k=1:height(meta)
    d=readtable(fullfile(root,meta.segment_csv(k)));
    if abs(d.angle_deg(1))>initial_limit
        error('Selected segment %s exceeds the retained display-offset rule.',meta.segment_id(k));
    end
    t=d.time_s-d.time_s(1); y=d.angle_deg-d.angle_deg(1);
    keep=t>=0 & t<=end_ms/1000; t=t(keep); y=y(keep);
    [t,ia]=unique(t,'stable'); y=y(ia);
    assert(numel(t)>=2);
    valid=dense_s>=t(1) & dense_s<=t(end);
    yvalid=interp1(t,y,dense_s(valid),'pchip');
    win=min(30,numel(yvalid)); if mod(win,2)==0, win=win-1; end
    yvalid=smoothdata(yvalid,'gaussian',win);
    yvalid=yvalid-yvalid(1);
    ydense=nan(size(dense_s)); ydense(valid)=yvalid;
    ydense(dense_s>t(end))=yvalid(end); ydense(dense_s<t(1))=yvalid(1);
    display.(char(meta.segment_id(k)))=ydense;
    plot(ax,dense_ms,ydense,'LineWidth',2,'Color',colors(k,:));
    text(ax,end_ms+2,ydense(end),strrep(meta.segment_id(k),'_',' '), ...
        'Color',colors(k,:),'FontSize',11,'Interpreter','none');
end
xlabel(ax,'Time from first retained sample (ms)');
ylabel(ax,sprintf('%s displacement (deg)',char(axis_name)));
title(ax,sprintf('Selected %s trajectories: display reconstruction',char(axis_name)));
xlim(ax,[0 end_ms+30]); set(ax,'FontSize',12);
stem=fullfile(output_dir,"figure3_"+axis_name+"_repository_replot");
exportgraphics(fig,stem+".pdf",'ContentType','vector');
exportgraphics(fig,stem+".png",'Resolution',180);
writetable(display,stem+"_display_samples.csv"); close(fig);
end

function plot_response_poles(fits,axis_name,output_dir)
fig=figure('Color','w','Position',[100 100 850 620]);
ax=axes(fig); hold(ax,'on'); grid(ax,'on'); colors=lines(height(fits));
handles=gobjects(height(fits),1);
for k=1:height(fits)
    xr=fits.pole_real_per_s(k); yi=fits.pole_imag_abs_per_s(k);
    handles(k)=plot(ax,[xr xr],[-yi yi],'o-','Color',colors(k,:),'LineWidth',1.4);
end
xmin=1.20*min(fits.pole_real_per_s); xx=linspace(xmin,0,400);
for zeta=[.6 .9]
    yy=-xx*sqrt(1-zeta^2)/zeta;
    plot(ax,xx,yy,'--','Color',[.35 .35 .35]);
    plot(ax,xx,-yy,'--','Color',[.35 .35 .35]);
end
xline(ax,0,'k-'); yline(ax,0,'k-'); axis(ax,'equal');
xlim(ax,[xmin 5]); ymax=1.2*max(fits.pole_imag_abs_per_s); ylim(ax,[-ymax ymax]);
legend(ax,handles,strrep(string(fits.segment_id),'_',' '),'Location','eastoutside','Interpreter','none');
xlabel(ax,'Real part (s^{-1})'); ylabel(ax,'Imaginary part (s^{-1})');
title(ax,sprintf('%s poles: repository refit; damping guides 0.6 and 0.9',char(axis_name)));
set(ax,'FontSize',12);
stem=fullfile(output_dir,"figure4_"+axis_name+"_repository_refit");
exportgraphics(fig,stem+".pdf",'ContentType','vector');
exportgraphics(fig,stem+".png",'Resolution',180); close(fig);
end

function write_json(path,value)
fid=fopen(path,'w'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end
