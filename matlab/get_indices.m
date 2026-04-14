function index = get_indices(n_bus,n_macs,n_branches,n_shunts,opt,n_gfl)
% usage: index = get_indices(n_bus,n_macs,n_branches,n_shunts,opt,n_gfl)
% produces a structure that helps us to know where stuff is in x,y,f,g
if nargin < 6
    n_gfl = 0;
end
if nargin < 5
    opt = psoptions;
end

interleave = opt.sim.interleave;    % default = 'true' or '1'

% differential (generator,exciter,governor,relays) variable index
nxsmac = 7;				% number of differential variables per machine 
                        % temperature is a fake differential variable. 
if interleave
    index.x.delta       = (1:nxsmac:nxsmac*n_macs);
    index.x.omega_pu    = (2:nxsmac:nxsmac*n_macs);
    index.x.Pm          = (3:nxsmac:nxsmac*n_macs);
    index.x.Eap         = (4:nxsmac:nxsmac*n_macs);
	index.x.E1 			= (5:nxsmac:nxsmac*n_macs);
	index.x.Efd			= (6:nxsmac:nxsmac*n_macs);
    index.x.P3			= (7:nxsmac:nxsmac*n_macs);
    % --- GFL states (new) ---
    gfl_start = nxsmac*n_macs;
    index.x.rho_gfl     = (1:n_gfl) + gfl_start;
    index.x.xi_pll      = (1:n_gfl) + gfl_start + n_gfl;
    index.x.xi_id       = (1:n_gfl) + gfl_start + 2*n_gfl;
    index.x.xi_iq       = (1:n_gfl) + gfl_start + 3*n_gfl;
    index.x.id_gfl      = (1:n_gfl) + gfl_start + 4*n_gfl;
    index.x.iq_gfl      = (1:n_gfl) + gfl_start + 5*n_gfl;
    % keep branch temperature at the end
    index.x.temp        = (1:n_branches) + gfl_start + 6*n_gfl;
else
    index.x.delta       = (1:n_macs);
    index.x.omega_pu    = (1:n_macs) + n_macs;
    index.x.Pm          = (1:n_macs) + n_macs*2;
    index.x.Eap         = (1:n_macs) + n_macs*3; 
	index.x.E1			= (1:n_macs) + n_macs*4;
	index.x.Efd 		= (1:n_macs) + n_macs*5; 
    index.x.P3  		= (1:n_macs) + n_macs*6;
    % --- GFL states (new) ---
    gfl_start = n_macs*7;
    index.x.rho_gfl     = (1:n_gfl) + gfl_start;
    index.x.xi_pll      = (1:n_gfl) + gfl_start + n_gfl;
    index.x.xi_id       = (1:n_gfl) + gfl_start + 2*n_gfl;
    index.x.xi_iq       = (1:n_gfl) + gfl_start + 3*n_gfl;
    index.x.id_gfl      = (1:n_gfl) + gfl_start + 4*n_gfl;
    index.x.iq_gfl      = (1:n_gfl) + gfl_start + 5*n_gfl;
    % keep branch temperature at the end
    index.x.temp        = (1:n_branches) + gfl_start + 6*n_gfl;
end

index.nx            = nxsmac*n_macs + 6*n_gfl + n_branches;
% index.nx            = nxsmac*n_macs;
index.x.omega       = index.x.omega_pu; 

% differential equation index is the same as x index
index.f.delta_dot = index.x.delta;
index.f.omega_dot = index.x.omega_pu;
index.f.Pm_dot    = index.x.Pm;
index.f.Eap_dot   = index.x.Eap;
index.f.E1_dot 	  = index.x.E1;
index.f.Efd_dot   = index.x.Efd;
index.f.P3_dot    = index.x.P3;
% --- GFL differential equations (new) ---
index.f.rho_gfl_dot = index.x.rho_gfl;
index.f.xi_pll_dot  = index.x.xi_pll;
index.f.xi_id_dot   = index.x.xi_id;
index.f.xi_iq_dot   = index.x.xi_iq;
index.f.id_gfl_dot  = index.x.id_gfl;
index.f.iq_gfl_dot  = index.x.iq_gfl;
% index.f.temp_dot  = index.x.temp;
% index.nf = index.nx;
index.nf = index.nx - n_branches;   % exclude the fake differential variable, temperature
index.f.swing = index.f.omega_dot;
angle_ref = opt.sim.angle_ref;    % 0:delta_sys  1: delta_coi

if angle_ref
    % Using 'meshgrid' to locate the cross derivatives terms when using COI
    [OmeDot,Del] = meshgrid(index.f.omega_dot, index.x.delta);
    index.COI.omega_dot = OmeDot(:);
    index.COI.delta = Del(:);
    [DelDot,~] = meshgrid(index.f.delta_dot, index.x.omega_pu);
    index.COI.delta_dot = DelDot(:);
    [~,Omega] = meshgrid(index.f.omega_dot, index.x.omega_pu);
    index.COI.omega_pu = Omega(:);
    [Eap,~] = meshgrid(index.f.Eap_dot, index.x.delta);
    index.COI.Eap_dot = Eap(:);
end

% algebraic variable index 
if interleave
    index.y.Vmag    = (1:2:2*n_bus);
    index.y.theta   = (2:2:2*n_bus); 
else
    index.y.Vmag    = (1:n_bus);
    index.y.theta   = (1:n_bus) + n_bus;
end

if ~angle_ref
    index.y.delta_sys = n_bus*2 + 1;
    index.ny = 2*n_bus + 1;
else
    index.ny = 2*n_bus;
end

% algebraic function index
index.g.P       = index.y.Vmag;
index.g.Q       = index.y.theta;
if ~angle_ref
    index.g.slack   = index.y.delta_sys;
end

index.ng = index.ny;

% relay index
index.re.temp   = 1:n_branches;
index.re.oc = (1:n_branches) + n_branches;
index.re.uvls = (1:n_shunts) + n_branches + n_branches;
index.re.ufls = (1:n_shunts)  + n_branches + n_branches + n_shunts;
index.re.dist = (1:n_branches) + n_branches + n_branches + n_shunts + n_shunts;
% --- GFL relay blocks (new) ---
index.re.gfl_uv = (1:n_gfl) + 3*n_branches + 2*n_shunts;
index.re.gfl_oc = (1:n_gfl) + 3*n_branches + 2*n_shunts + n_gfl;
index.re.gfl_pll = (1:n_gfl) + 3*n_branches + 2*n_shunts + 2*n_gfl;
index.re.nrelay = 3*n_branches + 2*n_shunts + 3*n_gfl;
return
