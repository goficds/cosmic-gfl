function gfl_results = plot_gfl_results(outputs,ps,opt,varargin)
% usage: gfl_results = plot_gfl_results(outputs,ps,opt)
% minimal single-/multi-GFL result extraction and plotting helper

C = psconstants;

if nargin < 3 || isempty(opt)
    opt = psoptions;
end
if nargin < 2 || isempty(ps)
    error('plot_gfl_results:err','ps is required.');
end
if nargin < 1 || isempty(outputs) || ~isfield(outputs,'outfilename') || isempty(outputs.outfilename)
    error('plot_gfl_results:err','outputs.outfilename is required.');
end
if ~isfield(ps,'gfl') || isempty(ps.gfl)
    error('plot_gfl_results:err','ps.gfl is required.');
end

% options
selected_gfl = 1;
save_plots = true;
visible = 'off';
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'gfl_id'
            selected_gfl = varargin{k+1};
        case 'save_plots'
            save_plots = varargin{k+1};
        case 'visible'
            visible = varargin{k+1};
        otherwise
            error('plot_gfl_results:err','unknown option: %s',varargin{k});
    end
end

% choose GFL row by idnum when possible
if isscalar(selected_gfl) && any(ps.gfl(:,C.gfl.idnum) == selected_gfl)
    gfl_row = find(ps.gfl(:,C.gfl.idnum) == selected_gfl,1,'first');
else
    gfl_row = selected_gfl;
end
if gfl_row < 1 || gfl_row > size(ps.gfl,1)
    error('plot_gfl_results:err','requested GFL index/id is out of range.');
end

data = readmatrix(outputs.outfilename,'OutputType','double','FileType','text','NumHeaderLines',1);
if isempty(data)
    error('plot_gfl_results:err','output file %s is empty.',outputs.outfilename);
end

n  = size(ps.bus,1);
ng = size(ps.mac,1);
m  = size(ps.branch,1);
n_sh = size(ps.shunt,1);
n_gfl = size(ps.gfl,1);
ix = get_indices(n,ng,m,n_sh,opt,n_gfl);

if size(data,2) < (1 + ix.nx + ix.ny)
    error('plot_gfl_results:err','output file does not contain the expected state columns.');
end

t = data(:,1);
X = data(:,1+(1:ix.nx));
Y = data(:,1+ix.nx+(1:ix.ny));

bus_no = ps.gfl(gfl_row,C.gfl.bus);
bus_i = ps.bus_i(bus_no);

rho = X(:,ix.x.rho_gfl(gfl_row));
xi_pll = X(:,ix.x.xi_pll(gfl_row));
xi_id = X(:,ix.x.xi_id(gfl_row));
xi_iq = X(:,ix.x.xi_iq(gfl_row));
id = X(:,ix.x.id_gfl(gfl_row));
iq = X(:,ix.x.iq_gfl(gfl_row));
vmag = Y(:,ix.y.Vmag(bus_i));
theta = Y(:,ix.y.theta(bus_i));

% event extraction
trip_close_mask = [];
event_times = [];
event_labels = {};
if isfield(outputs,'event_record') && ~isempty(outputs.event_record)
    event_times = outputs.event_record(:,C.ev.time);
    event_labels = cell(size(event_times));
    for ii = 1:numel(event_times)
        event_labels{ii} = local_event_label(outputs.event_record(ii,:),C);
    end
    trip_close_mask = outputs.event_record(:,C.ev.type) == C.ev.trip_branch | ...
                      outputs.event_record(:,C.ev.type) == C.ev.close_branch | ...
                      outputs.event_record(:,C.ev.type) == C.ev.trip_bus | ...
                      outputs.event_record(:,C.ev.type) == C.ev.trip_gfl | ...
                      outputs.event_record(:,C.ev.type) == C.ev.gfl_set_pref;
else
    trip_close_mask = false(0,1);
end

fig = figure('Name',sprintf('GFL %d results',selected_gfl), ...
    'Color','w','Visible',visible,'Position',[100 100 1100 850]);

ax1 = subplot(3,2,1);
plot(t,vmag,'LineWidth',1.2);
ylabel('Vmag (pu)');
title(sprintf('Bus %d voltage',bus_no));
grid on;
local_mark_events(ax1,event_times,event_labels,trip_close_mask);

ax2 = subplot(3,2,2);
plot(t,theta*180/pi,'LineWidth',1.2);
ylabel('\theta (deg)');
title(sprintf('Bus %d angle',bus_no));
grid on;
local_mark_events(ax2,event_times,event_labels,trip_close_mask);

ax3 = subplot(3,2,3);
plot(t,rho,'LineWidth',1.2);
ylabel('\rho (rad)');
title('PLL angle state');
grid on;
local_mark_events(ax3,event_times,event_labels,trip_close_mask);

ax4 = subplot(3,2,4);
plot(t,[id iq],'LineWidth',1.2);
ylabel('Current (pu)');
title('Current states');
legend({'i_d','i_q'},'Location','best');
grid on;
local_mark_events(ax4,event_times,event_labels,trip_close_mask);

ax5 = subplot(3,2,5);
plot(t,xi_pll,'LineWidth',1.2);
ylabel('\xi_{pll}');
xlabel('Time (s)');
title('PLL integrator');
grid on;
local_mark_events(ax5,event_times,event_labels,trip_close_mask);

ax6 = subplot(3,2,6);
plot(t,[xi_id xi_iq],'LineWidth',1.2);
ylabel('Integrator state');
xlabel('Time (s)');
title('Current-loop integrators');
legend({'\xi_{id}','\xi_{iq}'},'Location','best');
grid on;
local_mark_events(ax6,event_times,event_labels,trip_close_mask);

sgtitle(sprintf('GFL %d at bus %d',selected_gfl,bus_no));

% package extracted results
[gfl_dir,gfl_name,~] = fileparts(outputs.outfilename);
if isempty(gfl_dir)
    gfl_dir = pwd;
end
png_file = fullfile(gfl_dir,[gfl_name '_gfl' num2str(selected_gfl) '.png']);
fig_file = fullfile(gfl_dir,[gfl_name '_gfl' num2str(selected_gfl) '.fig']);

gfl_results = struct();
gfl_results.t = t;
gfl_results.rho = rho;
gfl_results.xi_pll = xi_pll;
gfl_results.xi_id = xi_id;
gfl_results.xi_iq = xi_iq;
gfl_results.id = id;
gfl_results.iq = iq;
gfl_results.Vmag = vmag;
gfl_results.theta = theta;
gfl_results.bus_no = bus_no;
gfl_results.gfl_row = gfl_row;
gfl_results.event_times = event_times;
gfl_results.event_labels = event_labels;
gfl_results.figure_handle = fig;
gfl_results.png_file = png_file;
gfl_results.fig_file = fig_file;

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
