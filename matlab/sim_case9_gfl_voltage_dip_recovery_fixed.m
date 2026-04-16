%% case9: single GFL with fixed branch trip/close voltage dip-recovery proxy
clearvars; clear global; close all; clc; C = psconstants;

% do not touch path if we are deploying code
if ~(ismcc || isdeployed)
    addpath('../data');
    addpath('../numerics');
end

% simulation time
trip_t = 0.2;
close_t = 0.5;
t_max = 1.5;
branch_id = 6;      % branch 5-7
monitor_bus = 5;    % track the nearby load bus voltage

% select data case
ps = updateps(case9_ps);
ps.branch(:,C.br.tap)       = 1;
ps.shunt(:,C.sh.frac_S)     = 1;
ps.shunt(:,C.sh.frac_E)     = 0;
ps.shunt(:,C.sh.frac_Z)     = 0;
ps.shunt(:,C.sh.gamma)      = 0.08;

% options
opt = psoptions;
opt.sim.integration_scheme = 1;
opt.sim.dt_default = 1e-3;
% discrete branch switching is more robust with Newton warm-start than fsolve here
opt.nr.use_fsolve = false;
opt.verbose = true;
opt.sim.gen_control = 1;
opt.sim.angle_ref = 0;
opt.sim.COI_weight = 0;
opt.sim.trip_gfl_only_island = false;

% initialize solved operating point
ps = newpf(ps,opt);
[ps.Ybus,ps.Yf,ps.Yt] = getYbus(ps,false);
ps = update_load_freq_source(ps);
[ps.mac,ps.exc,ps.gov] = get_mac_state(ps,'salient');

% add one online GFL (small active-power reference)
gfl = zeros(1,C.gfl.cols);
gfl_bus_no = ps.bus(1,C.bu.id);
gfl_bus_i = ps.bus_i(gfl_bus_no);
gfl(C.gfl.bus) = gfl_bus_no;
gfl(C.gfl.status) = 1;
gfl(C.gfl.Sn) = 100;
gfl(C.gfl.Pref) = 1;
gfl(C.gfl.Qref) = 0;
gfl(C.gfl.Kp_pll) = 10;
gfl(C.gfl.Ki_pll) = 200;
gfl(C.gfl.Kp_i) = 2;
gfl(C.gfl.Ki_i) = 50;
gfl(C.gfl.Lf) = 0.1;
gfl(C.gfl.Rf) = 0.01;
gfl(C.gfl.Imax) = 2;
gfl(C.gfl.Vtrip) = 0.5;
gfl(C.gfl.wmax_dev) = 2*pi*5;
% initialize after power-flow using solved bus angle
gfl(C.gfl.rho) = ps.bus(gfl_bus_i,C.bu.Vang)*pi/180;
gfl(C.gfl.xi_pll) = 0;
gfl(C.gfl.xi_id) = 0;
gfl(C.gfl.xi_iq) = 0;
gfl(C.gfl.id) = 0;
gfl(C.gfl.iq) = 0;
gfl(C.gfl.idnum) = 1;
ps.gfl = gfl;
ps = updateps(ps);

% initialize relays
ps.relay = get_relays(ps,'all',opt);

% initialize global relay states
global t_delay t_prev_check dist2threshold state_a gfl_rocof_state
n    = size(ps.bus,1);
ng   = size(ps.mac,1);
m    = size(ps.branch,1);
n_sh = size(ps.shunt,1);
n_gfl = size(ps.gfl,1);
ix   = get_indices(n,ng,m,n_sh,opt,n_gfl);
t_delay = inf(size(ps.relay,1),1);
t_delay([ix.re.uvls])= opt.sim.uvls_tdelay_ini;
t_delay([ix.re.ufls])= opt.sim.ufls_tdelay_ini;
t_delay([ix.re.dist])= opt.sim.dist_tdelay_ini;
t_delay([ix.re.temp])= opt.sim.temp_tdelay_ini;
t_prev_check = nan(size(ps.relay,1),1);
dist2threshold = inf(size(ix.re.oc,2)*2,1);
state_a = zeros(size(ix.re.oc,2)*2,1);

% start + fixed branch trip + fixed branch close + finish
event = zeros(4,C.ev.cols);
event(1,[C.ev.time C.ev.type]) = [0 C.ev.start];
event(2,[C.ev.time C.ev.type C.ev.branch_loc]) = [trip_t C.ev.trip_branch branch_id];
event(3,[C.ev.time C.ev.type C.ev.branch_loc]) = [close_t C.ev.close_branch branch_id];
event(4,[C.ev.time C.ev.type]) = [t_max C.ev.finish];

% run simulation
[outputs,ps_end] = simgrid(ps,event,'sim_case9_gfl_voltage_dip_recovery_fixed',opt);

% completion check: use outputs.t_simulated because CSV tails may be sparse after switching
if ~outputs.success || ~isfield(outputs,'t_simulated') || isempty(outputs.t_simulated) || outputs.t_simulated(end) < (t_max - 1e-6)
    error('sim_case9_gfl_voltage_dip_recovery_fixed failed: t_simulated_end = %.6f, success = %d',outputs.t_simulated(end),outputs.success);
end

% verify that both discrete actions were recorded
if ~isfield(outputs,'event_record') || isempty(outputs.event_record)
    error('sim_case9_gfl_voltage_dip_recovery_fixed failed: no events recorded.');
end
has_trip = any(outputs.event_record(:,C.ev.type) == C.ev.trip_branch & outputs.event_record(:,C.ev.branch_loc) == branch_id);
has_close = any(outputs.event_record(:,C.ev.type) == C.ev.close_branch & outputs.event_record(:,C.ev.branch_loc) == branch_id);
if ~has_trip || ~has_close
    error('sim_case9_gfl_voltage_dip_recovery_fixed failed: missing trip/close event record for branch %d.',branch_id);
end

% inspect voltage dip/recovery using algebraic proxy states around the switching events
[x0,y0] = get_xy(ps,opt);
v_pre = y0(ix.y.Vmag(ps.bus_i(monitor_bus)));

ps_trip = ps;
ps_trip.event_record = [];
[ps_trip,~] = process_event(ps_trip,event(2,:),opt);
[ps_trip.Ybus,ps_trip.Yf,ps_trip.Yt,ps_trip.Yft,ps_trip.sh_ft] = getYbus(ps_trip,false);
y_trip = solve_algebraic(trip_t,x0,y0,ps_trip,opt);
if isempty(y_trip)
    error('sim_case9_gfl_voltage_dip_recovery_fixed failed: could not solve algebraic state after branch trip.');
end
v_dip = y_trip(ix.y.Vmag(ps_trip.bus_i(monitor_bus)));

ps_close = ps_trip;
ps_close.event_record = [];
[ps_close,~] = process_event(ps_close,event(3,:),opt);
[ps_close.Ybus,ps_close.Yf,ps_close.Yt,ps_close.Yft,ps_close.sh_ft] = getYbus(ps_close,false);
y_close = solve_algebraic(close_t,x0,y_trip,ps_close,opt);
if isempty(y_close)
    error('sim_case9_gfl_voltage_dip_recovery_fixed failed: could not solve algebraic state after branch close.');
end
v_rec = y_close(ix.y.Vmag(ps_close.bus_i(monitor_bus)));

if ~(v_dip < v_pre - 1e-4)
    error('sim_case9_gfl_voltage_dip_recovery_fixed failed: no voltage dip observed at bus %d (pre = %.6f, dip = %.6f).',monitor_bus,v_pre,v_dip);
end
if ~(abs(v_rec - v_pre) <= max(5e-3,0.5*(v_pre - v_dip)))
    error('sim_case9_gfl_voltage_dip_recovery_fixed failed: voltage did not recover at bus %d (pre = %.6f, rec = %.6f).',monitor_bus,v_pre,v_rec);
end

gfl_results = plot_gfl_results(outputs,ps_end,opt,'gfl_id',1,'visible','off','save_plots',true);
if ~isfile(gfl_results.png_file)
    error('sim_case9_gfl_voltage_dip_recovery_fixed failed: GFL plot image was not created.');
end
fprintf(['sim_case9_gfl_voltage_dip_recovery_fixed passed: reached t = %.6f s, ' ...
    'recorded branch %d trip/close, Vbus%d %.6f -> %.6f -> %.6f pu; plot saved to %s.\n'], ...
    outputs.t_simulated(end),branch_id,monitor_bus,v_pre,v_dip,v_rec,gfl_results.png_file);
