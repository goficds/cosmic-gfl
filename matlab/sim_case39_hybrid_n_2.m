function [outputs,ps_end,plot_results] = sim_case39_hybrid_n_2(a,b)
%% simulate 39-bus hybrid case with synchronous machines and one GFL
if nargin < 1 || isempty(a)
    a = 35;
end
if nargin < 2 || isempty(b)
    b = 23;
end

clearvars -except a b; clear global; close all; clc; C = psconstants;

if ~(ismcc || isdeployed)
    addpath('../data');
    addpath('../numerics');
end

t_max = 40;

ps = updateps(case39_ps);
ps = unify_generators(ps);
ps.branch(:,C.br.tap)       = 1;
ps.shunt(:,C.sh.factor)     = 1;
ps.shunt(:,C.sh.status)     = 1;
ps.shunt(:,C.sh.frac_S)     = 1;
ps.shunt(:,C.sh.frac_E)     = 0;
ps.shunt(:,C.sh.frac_Z)     = 0;
ps.shunt(:,C.sh.gamma)      = 0.08;

rateB_rateA = ps.branch(:,C.br.rateB)./ps.branch(:,C.br.rateA);
rateC_rateA = ps.branch(:,C.br.rateC)./ps.branch(:,C.br.rateA);
ps.branch(rateB_rateA==1,C.br.rateB) = 1.1 * ps.branch(rateB_rateA==1,C.br.rateA);
ps.branch(rateC_rateA==1,C.br.rateC) = 1.5 * ps.branch(rateC_rateA==1,C.br.rateA);

opt = psoptions;
opt.sim.integration_scheme = 1;
opt.sim.dt_default = 1/20;
opt.nr.use_fsolve = false;
opt.verbose = true;
opt.sim.gen_control = 1;
opt.sim.angle_ref = 0;
opt.sim.COI_weight = 0;
opt.sim.uvls_tdelay_ini = 0.5;
opt.sim.ufls_tdelay_ini = 0.5;
opt.sim.dist_tdelay_ini = 0.5;
opt.sim.temp_tdelay_ini = 0;
opt.sim.trip_gfl_only_island = false;

ps = newpf(ps,opt);
[ps.Ybus,ps.Yf,ps.Yt] = getYbus(ps,false);
ps = update_load_freq_source(ps);
[ps.mac,ps.exc,ps.gov] = get_mac_state(ps,'salient');

% add one hybrid GFL resource at a load bus
bus_no = 4;
bus_i = ps.bus_i(bus_no);
gfl = zeros(1,C.gfl.cols);
gfl(C.gfl.bus) = bus_no;
gfl(C.gfl.status) = 1;
gfl(C.gfl.Sn) = 150;
gfl(C.gfl.Pref) = 40;
gfl(C.gfl.Qref) = 0;
gfl(C.gfl.Kp_pll) = 10;
gfl(C.gfl.Ki_pll) = 200;
gfl(C.gfl.Kp_i) = 2;
gfl(C.gfl.Ki_i) = 50;
gfl(C.gfl.Lf) = 0.1;
gfl(C.gfl.Rf) = 0.01;
gfl(C.gfl.Imax) = 2;
gfl(C.gfl.Vtrip) = 0.55;
gfl(C.gfl.wmax_dev) = 2*pi*5;
gfl(C.gfl.rho) = ps.bus(bus_i,C.bu.Vang)*pi/180;
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
t_delay(ix.re.uvls)= opt.sim.uvls_tdelay_ini;
t_delay(ix.re.ufls)= opt.sim.ufls_tdelay_ini;
t_delay(ix.re.dist)= opt.sim.dist_tdelay_ini;
t_delay(ix.re.temp)= opt.sim.temp_tdelay_ini;
t_delay(ix.re.gfl_rocof) = opt.sim.gfl_rocof_tdelay_ini;
t_prev_check = nan(size(ps.relay,1),1);
dist2threshold = inf(max(size(ps.relay,1),size(ix.re.oc,2)*2),1);
state_a = zeros(max(size(ps.relay,1),size(ix.re.oc,2)*2),1);
gfl_rocof_state = [];

event = zeros(4,C.ev.cols);
event(1,[C.ev.time C.ev.type]) = [0 C.ev.start];
event(2,[C.ev.time C.ev.type C.ev.branch_loc]) = [5 C.ev.trip_branch a];
event(3,[C.ev.time C.ev.type C.ev.branch_loc]) = [5 C.ev.trip_branch b];
event(4,[C.ev.time C.ev.type]) = [t_max C.ev.finish];

[outputs,ps_end] = simgrid(ps,event,'sim_case39_hybrid_n_2',opt);

if ~outputs.success || ~isfield(outputs,'t_simulated') || isempty(outputs.t_simulated)
    error('sim_case39_hybrid_n_2 failed: missing simulation completion flag or time vector.');
end
plot_results = plot_case39_hybrid_results(outputs,ps_end,opt,'visible','off','save_plots',true,'gfl_id',1);
if ~isfile(plot_results.png_file)
    error('sim_case39_hybrid_n_2 failed: hybrid plot image was not created.');
end
fprintf('sim_case39_hybrid_n_2 passed: reached t = %.6f s; hybrid plot saved to %s\n',outputs.t_simulated(end),plot_results.png_file);
end
