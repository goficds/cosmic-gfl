%% case9: single GFL with forced disturbance and ROCOF trip validation
clearvars; clear global; close all; clc; C = psconstants;

if ~(ismcc || isdeployed)
    addpath('../data');
    addpath('../numerics');
end

t_max = 1;
pref_step_t = 0.2;
branch_trip_t = 0.21;
branch_id = 6;

ps = updateps(case9_ps);
ps.branch(:,C.br.tap)       = 1;
ps.shunt(:,C.sh.frac_S)     = 1;
ps.shunt(:,C.sh.frac_E)     = 0;
ps.shunt(:,C.sh.frac_Z)     = 0;
ps.shunt(:,C.sh.gamma)      = 0.08;

opt = psoptions;
opt.sim.integration_scheme = 1;
opt.sim.dt_default = 1e-3;
opt.nr.use_fsolve = false;
opt.verbose = true;
opt.sim.gen_control = 1;
opt.sim.angle_ref = 0;
opt.sim.COI_weight = 0;
opt.sim.trip_gfl_only_island = false;
opt.sim.gfl_rocof_limit = 1e-6;
opt.sim.gfl_rocof_tdelay_ini = 0.0;

ps = newpf(ps,opt);
[ps.Ybus,ps.Yf,ps.Yt] = getYbus(ps,false);
ps = update_load_freq_source(ps);
[ps.mac,ps.exc,ps.gov] = get_mac_state(ps,'salient');

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
gfl(C.gfl.rho) = ps.bus(gfl_bus_i,C.bu.Vang)*pi/180;
gfl(C.gfl.xi_pll) = 0;
gfl(C.gfl.xi_id) = 0;
gfl(C.gfl.xi_iq) = 0;
gfl(C.gfl.id) = 0;
gfl(C.gfl.iq) = 0;
gfl(C.gfl.idnum) = 1;
ps.gfl = gfl;
ps = updateps(ps);

ps.relay = get_relays(ps,'all',opt);

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
t_delay([ix.re.gfl_rocof]) = opt.sim.gfl_rocof_tdelay_ini;
t_prev_check = nan(size(ps.relay,1),1);
dist2threshold = inf(max(size(ps.relay,1),size(ix.re.oc,2)*2),1);
state_a = zeros(max(size(ps.relay,1),size(ix.re.oc,2)*2),1);
gfl_rocof_state = [];

event = zeros(4,C.ev.cols);
event(1,[C.ev.time C.ev.type]) = [0 C.ev.start];
event(2,[C.ev.time C.ev.type C.ev.quantity C.ev.gfl_loc]) = [pref_step_t C.ev.gfl_set_pref 10 1];
event(3,[C.ev.time C.ev.type C.ev.branch_loc]) = [branch_trip_t C.ev.trip_branch branch_id];
event(4,[C.ev.time C.ev.type]) = [t_max C.ev.finish];

[outputs,ps_end] = simgrid(ps,event,'sim_case9_gfl_rocof_fixed',opt);

if ~outputs.success || ~isfield(outputs,'t_simulated') || isempty(outputs.t_simulated) || outputs.t_simulated(end) < (t_max - 1e-6)
    error('sim_case9_gfl_rocof_fixed failed: t_simulated_end = %.6f, success = %d',outputs.t_simulated(end),outputs.success);
end
if ~isfield(outputs,'event_record') || isempty(outputs.event_record)
    error('sim_case9_gfl_rocof_fixed failed: no events recorded.');
end
has_pref_step = any(outputs.event_record(:,C.ev.type) == C.ev.gfl_set_pref & outputs.event_record(:,C.ev.gfl_loc) == 1);
has_branch_trip = any(outputs.event_record(:,C.ev.type) == C.ev.trip_branch & outputs.event_record(:,C.ev.branch_loc) == branch_id);
has_gfl_trip = any(outputs.event_record(:,C.ev.type) == C.ev.trip_gfl & outputs.event_record(:,C.ev.gfl_loc) == 1);
if ~has_pref_step
    error('sim_case9_gfl_rocof_fixed failed: GFL Pref step was not recorded.');
end
if ~has_branch_trip
    error('sim_case9_gfl_rocof_fixed failed: fixed branch trip was not recorded.');
end
if ~has_gfl_trip
    error('sim_case9_gfl_rocof_fixed failed: ROCOF-driven GFL trip was not recorded.');
end

gfl_results = plot_gfl_results(outputs,ps_end,opt,'gfl_id',1,'visible','off','save_plots',true);
if ~isfile(gfl_results.png_file)
    error('sim_case9_gfl_rocof_fixed failed: GFL plot image was not created.');
end
fprintf('sim_case9_gfl_rocof_fixed passed: reached t = %.6f s with forced disturbance and recorded GFL ROCOF trip; plot saved to %s\n',outputs.t_simulated(end),gfl_results.png_file);
