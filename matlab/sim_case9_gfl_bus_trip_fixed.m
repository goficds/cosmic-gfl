%% case9: single GFL with fixed bus-trip disturbance
clearvars; clear global; close all; clc; C = psconstants;

% do not touch path if we are deploying code
if ~(ismcc || isdeployed)
    addpath('../data');
    addpath('../numerics');
end

% simulation time
t_max = 1;

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
% after a discrete event, Newton warm-start is typically more robust than fsolve
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
gfl_bus_i = ps.bus_i(gfl_bus_no); %#ok<NASGU>
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
gfl(C.gfl.rho) = ps.bus(ps.bus_i(gfl_bus_no),C.bu.Vang)*pi/180;
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
global t_delay t_prev_check dist2threshold state_a
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

% start + fixed bus trip + finish (single deterministic discrete action)
event = zeros(3,C.ev.cols);
event(1,[C.ev.time C.ev.type]) = [0 C.ev.start];
event(2,[C.ev.time C.ev.type]) = [0.2 C.ev.trip_bus];
% choose a load-only bus to avoid creating a generator-only island
event(2,C.ev.bus_loc) = 5;
event(3,[C.ev.time C.ev.type]) = [t_max C.ev.finish];

% run simulation
[outputs,~] = simgrid(ps,event,'sim_case9_gfl_bus_trip_fixed',opt);

% dedicated smoke-test post path
% Note: for partitioned networks, the CSV tail can contain NaNs and stop updating.
% Use outputs.t_simulated as the primary completion signal.
t_last = NaN;
if isfield(outputs,'outfilename') && ~isempty(outputs.outfilename)
    data = readmatrix(outputs.outfilename,'OutputType','double','FileType','text','NumHeaderLines',1);
    if ~isempty(data)
        t_last = data(end,1);
    end
end
if ~outputs.success || ~isfield(outputs,'t_simulated') || isempty(outputs.t_simulated) || outputs.t_simulated(end) < (t_max - 1e-6)
    error('sim_case9_gfl_bus_trip_fixed failed: t_simulated_end = %.6f, success = %d',outputs.t_simulated(end),outputs.success);
end

% verify that the bus-trip discrete action was recorded
if ~isfield(outputs,'event_record') || isempty(outputs.event_record)
    error('sim_case9_gfl_bus_trip_fixed failed: no events recorded.');
end
if ~any(outputs.event_record(:,C.ev.type) == C.ev.trip_bus)
    error('sim_case9_gfl_bus_trip_fixed failed: fixed bus-trip action was not recorded.');
end

fprintf('sim_case9_gfl_bus_trip_fixed passed: reached t = %.6f s with recorded bus trip event.\n',outputs.t_simulated(end));
