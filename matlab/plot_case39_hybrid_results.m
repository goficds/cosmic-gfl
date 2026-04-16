function hybrid_results = plot_case39_hybrid_results(outputs,ps,opt,varargin)
% usage: hybrid_results = plot_case39_hybrid_results(outputs,ps,opt,varargin)
% combined plotting for synchronous-machine and GFL continuous states

C = psconstants;

if nargin < 3 || isempty(opt)
    opt = psoptions;
end
if nargin < 2 || isempty(ps)
    error('plot_case39_hybrid_results:err','ps is required.');
end
if nargin < 1 || isempty(outputs) || ~isfield(outputs,'outfilename') || isempty(outputs.outfilename)
    error('plot_case39_hybrid_results:err','outputs.outfilename is required.');
end

visible = 'off';
save_plots = true;
selected_gfl = 1;
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'visible'
            visible = varargin{k+1};
        case 'save_plots'
            save_plots = varargin{k+1};
        case 'gfl_id'
            selected_gfl = varargin{k+1};
        otherwise
            error('plot_case39_hybrid_results:err','unknown option: %s',varargin{k});
    end
end

data = readmatrix(outputs.outfilename,'OutputType','double','FileType','text','NumHeaderLines',1);
if isempty(data)
    error('plot_case39_hybrid_results:err','output file %s is empty.',outputs.outfilename);
end

n  = size(ps.bus,1);
ng = size(ps.mac,1);
m  = size(ps.branch,1);
n_sh = size(ps.shunt,1);
n_gfl = 0;
if isfield(ps,'gfl') && ~isempty(ps.gfl)
    n_gfl = size(ps.gfl,1);
end
ix = get_indices(n,ng,m,n_sh,opt,n_gfl);
if size(data,2) < (1 + ix.nx + ix.ny)
    error('plot_case39_hybrid_results:err','output file does not contain the expected state columns.');
end

t = data(:,1);
X = data(:,1+(1:ix.nx));
Y = data(:,1+ix.nx+(1:ix.ny));

omega = X(:,ix.x.omega_pu);
delta = X(:,ix.x.delta);
Pm = X(:,ix.x.Pm);
Vmag = Y(:,ix.y.Vmag);

has_gfl = n_gfl > 0;
if has_gfl
    if isscalar(selected_gfl) && any(ps.gfl(:,C.gfl.idnum) == selected_gfl)
        gfl_row = find(ps.gfl(:,C.gfl.idnum) == selected_gfl,1,'first');
    else
        gfl_row = selected_gfl;
    end
    if gfl_row < 1 || gfl_row > n_gfl
        error('plot_case39_hybrid_results:err','requested GFL index/id is out of range.');
    end
    bus_no = ps.gfl(gfl_row,C.gfl.bus);
    bus_i = ps.bus_i(bus_no);
    rho = X(:,ix.x.rho_gfl(gfl_row));
    xi_pll = X(:,ix.x.xi_pll(gfl_row));
    id = X(:,ix.x.id_gfl(gfl_row));
    iq = X(:,ix.x.iq_gfl(gfl_row));
    Vgfl = Y(:,ix.y.Vmag(bus_i));
    Tgfl = Y(:,ix.y.theta(bus_i))*180/pi;
else
    gfl_row = [];
    bus_no = [];
    rho = [];
    xi_pll = [];
    id = [];
    iq = [];
    Vgfl = [];
    Tgfl = [];
end

event_times = [];
event_labels = {};
mark_mask = false(0,1);
if isfield(outputs,'event_record') && ~isempty(outputs.event_record)
    event_times = outputs.event_record(:,C.ev.time);
    event_labels = cell(size(event_times));
    for ii = 1:numel(event_times)
        event_labels{ii} = local_event_label(outputs.event_record(ii,:),C);
    end
    mark_mask = true(size(event_times));
end

fig = figure('Name','Case39 hybrid dynamics','Color','w','Visible',visible,'Position',[80 80 1400 900]);

a1 = subplot(3,2,1);
plot(t,omega,'LineWidth',1.0);
ylabel('\omega (pu)');
title('Synchronous-machine speed states');
grid on;
local_mark_events(a1,event_times,event_labels,mark_mask);

a2 = subplot(3,2,2);
plot(t,delta,'LineWidth',1.0);
ylabel('\delta (rad)');
title('Synchronous-machine angle states');
grid on;
local_mark_events(a2,event_times,event_labels,mark_mask);

a3 = subplot(3,2,3);
plot(t,Pm,'LineWidth',1.0);
ylabel('P_m (pu)');
xlabel('Time (s)');
title('Governor/mechanical power states');
grid on;
local_mark_events(a3,event_times,event_labels,mark_mask);

a4 = subplot(3,2,4);
plot(t,Vmag,'LineWidth',0.9);
ylabel('|V| (pu)');
xlabel('Time (s)');
title('All bus voltages');
grid on;
local_mark_events(a4,event_times,event_labels,mark_mask);

if has_gfl
    a5 = subplot(3,2,5);
    plot(t,[rho xi_pll],'LineWidth',1.1);
    ylabel('State');
    xlabel('Time (s)');
    title(sprintf('GFL %d PLL states at bus %d',ps.gfl(gfl_row,C.gfl.idnum),bus_no));
    legend({'rho','xi_{pll}'},'Location','best');
    grid on;
    local_mark_events(a5,event_times,event_labels,mark_mask);

    a6 = subplot(3,2,6);
    plot(t,[id iq Vgfl Tgfl],'LineWidth',1.1);
    ylabel('State / output');
    xlabel('Time (s)');
    title(sprintf('GFL %d current and bus quantities',ps.gfl(gfl_row,C.gfl.idnum)));
    legend({'i_d','i_q','V_{bus}','theta_{bus}(deg)'},'Location','best');
    grid on;
    local_mark_events(a6,event_times,event_labels,mark_mask);
else
    a5 = subplot(3,2,5);
    axis(a5,'off');
    text(0.1,0.5,'No GFL units in this case.','Parent',a5,'FontSize',12);
    a6 = subplot(3,2,6);
    axis(a6,'off');
end

sgtitle('Case39 hybrid synchronous-machine and GFL continuous states');

[out_dir,out_name,~] = fileparts(outputs.outfilename);
if isempty(out_dir)
    out_dir = pwd;
end
png_file = fullfile(out_dir,[out_name '_hybrid.png']);
fig_file = fullfile(out_dir,[out_name '_hybrid.fig']);

hybrid_results = struct();
hybrid_results.t = t;
hybrid_results.omega = omega;
hybrid_results.delta = delta;
hybrid_results.Pm = Pm;
hybrid_results.Vmag = Vmag;
hybrid_results.figure_handle = fig;
hybrid_results.png_file = png_file;
hybrid_results.fig_file = fig_file;
hybrid_results.event_times = event_times;
hybrid_results.event_labels = event_labels;
if has_gfl
    hybrid_results.gfl_row = gfl_row;
    hybrid_results.gfl_bus_no = bus_no;
    hybrid_results.rho = rho;
    hybrid_results.xi_pll = xi_pll;
    hybrid_results.id = id;
    hybrid_results.iq = iq;
    hybrid_results.Vgfl = Vgfl;
    hybrid_results.Tgfl = Tgfl;
end

if save_plots
    saveas(fig,png_file);
    savefig(fig,fig_file);
end
end

function local_mark_events(ax,event_times,event_labels,mark_mask)
axes(ax); %#ok<LAXES>
hold(ax,'on');
if isempty(event_times)
    return
end
for jj = 1:numel(event_times)
    if ~isempty(mark_mask) && numel(mark_mask) >= jj && ~mark_mask(jj)
        continue
    end
    xl = xline(ax,event_times(jj),'--','Color',[0.4 0.4 0.4],'LineWidth',0.8);
    if ~isempty(event_labels) && numel(event_labels) >= jj
        xl.Label = event_labels{jj};
        xl.LabelVerticalAlignment = 'middle';
        xl.LabelOrientation = 'aligned';
        xl.FontSize = 8;
    end
end
end

function label = local_event_label(event_row,C)
switch event_row(C.ev.type)
    case C.ev.start
        label = 'start';
    case C.ev.finish
        label = 'finish';
    case C.ev.trip_branch
        label = sprintf('trip br %d',event_row(C.ev.branch_loc));
    case C.ev.close_branch
        label = sprintf('close br %d',event_row(C.ev.branch_loc));
    case C.ev.trip_bus
        label = sprintf('trip bus %d',event_row(C.ev.bus_loc));
    case C.ev.trip_gfl
        label = sprintf('trip gfl %d',event_row(C.ev.gfl_loc));
    case C.ev.gfl_set_pref
        label = sprintf('Pref->%.3g',event_row(C.ev.quantity));
    otherwise
        label = sprintf('event %d',event_row(C.ev.type));
end
end
