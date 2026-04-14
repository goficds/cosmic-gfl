function [f,df_dx,df_dy] = differential_eqs(t,x,y,ps,opt) %#ok<INUSL>
% usage: [f,df_dx,df_dy] = differential_eqs(t,x,y,ps,opt)
% differential equations that model the elements of the power system
%
% inputs:
%   t   -> time in seconds
%   x   -> [delta omega Pm Eap E1 Efd] for each machine
%   y   -> [Vmag Theta]
%   ps  -> power system structure
%   opt -> options structure
%
% outputs: 
%   f(1) -> ddelta_dt
%   f(2) -> domega_dt
%   f(3) -> dPm_dt
%   f(4) -> dEap_dt
%   f(5) -> dE1_dt
%   f(6) -> dEfd_dt
%   f(7) -> dP3_dt

% constants
C           = psconstants;
n           = size(ps.bus,1);
ng          = size(ps.mac,1);
m           = size(ps.branch,1);
n_sh        = size(ps.shunt,1);
n_gfl       = 0;
if isfield(ps,'gfl') && ~isempty(ps.gfl)
    n_gfl = size(ps.gfl,1);
end
ix          = get_indices(n,ng,m,n_sh,opt,n_gfl);
gc          = opt.sim.gen_control; 
angle_ref = opt.sim.angle_ref;                 % angle reference: 0:delta_sys,1:delta_coi
COI_weight = opt.sim.COI_weight;               % weight of center of inertia

% extract parameters from ps
Xds     = ps.mac(:,C.ma.Xd);
Xdps    = ps.mac(:,C.ma.Xdp);
Xqs     = ps.mac(:,C.ma.Xq);
Td0ps   = ps.mac(:,C.ma.Td0p);
Ds      = ps.mac(:,C.ma.D);
Ms      = ps.mac(:,C.ma.M);
omega_0 = 2*pi*ps.frequency;

if COI_weight
    weight      = Ms;
else
    weight     = ps.gen(:,C.ge.mBase);
end

% extract differential variables
deltas      = x(ix.x.delta);
omegas_pu   = x(ix.x.omega_pu);
Pms         = x(ix.x.Pm);
Eaps        = x(ix.x.Eap);
Efds        = x(ix.x.Efd);
E1s         = x(ix.x.E1);
P3s         = x(ix.x.P3);

% extract algebraic variables and fix the slack bus angle
mac_buses   = ps.mac(:,C.mac.gen);
mac_bus_i   = ps.bus_i(mac_buses);
Vmags       = y(ix.y.Vmag);
Thetas      = y(ix.y.theta);
mac_Vmags   = Vmags(mac_bus_i);
mac_Thetas  = Thetas(mac_bus_i);

% machine angles, relative to the bus angles:
if ~angle_ref
    delta_sys = y(ix.y.delta_sys);
    delta_m = deltas + delta_sys - mac_Thetas;
else
    delta_coi = sum(weight.*deltas)/sum(weight);
    delta_m = deltas - delta_coi - mac_Thetas;
    omegas_coi = sum(weight.*omegas_pu,1)/sum(weight,1);
end


% calculate Pe: Eq. 7.81 from Bergen & Vittal
Pes = (Eaps.*mac_Vmags./Xdps).*sin(delta_m) + mac_Vmags.^2./2.*(1./Xqs-1./Xdps).*sin(2*delta_m);

% initialize output
f = zeros(length(x),1);

% calculate swing equations
if ~angle_ref
    f(ix.f.delta_dot) = omega_0.*(omegas_pu-1);
    f(ix.f.omega_dot) = (Pms - Pes - Ds.*(omegas_pu-1))./Ms;
else
    f(ix.f.delta_dot) = omega_0.*(omegas_pu-omegas_coi);
    f(ix.f.omega_dot) = (Pms - Pes - Ds.*(omegas_pu-omegas_coi))./Ms;
end
f(ix.f.Eap_dot)   = -Eaps.*Xds./(Td0ps.*Xdps)+(Xds./Xdps-1).*mac_Vmags.*cos(delta_m)./Td0ps+Efds./Td0ps; % Eq. 7.75 from Bergen & Vittal   

% calculate governor and exciter equations
if gc
    [f(ix.f.Pm_dot),f(ix.f.P3_dot),df_dx_gov]               = governor_eqs_modified([Pms';P3s'],omegas_pu,ps);
    [f(ix.f.Efd_dot),f(ix.f.E1_dot),df_dx_exc,df_dy_exc]   	= exciter_eqs([Efds';E1s'],mac_Vmags,ps);
end

% --- GFL dynamics (new) ---
if n_gfl > 0
    gfl_on = ps.gfl(:,C.gfl.status) > 0;
    if any(gfl_on)
        eps_v = 1e-4; % TODO: expose as an option
        rho = x(ix.x.rho_gfl);
        xi_pll = x(ix.x.xi_pll);
        xi_id = x(ix.x.xi_id);
        xi_iq = x(ix.x.xi_iq);
        id = x(ix.x.id_gfl);
        iq = x(ix.x.iq_gfl);
        gfl_bus_i = ps.bus_i(ps.gfl(:,C.gfl.bus));
        Vgfl = Vmags(gfl_bus_i);
        Tgfl = Thetas(gfl_bus_i);
        delta_theta = Tgfl - rho;
        vd = Vgfl .* cos(delta_theta);
        vq = Vgfl .* sin(delta_theta);
        vd_eff = max(vd,eps_v);

        Pref = ps.gfl(:,C.gfl.Pref) ./ max(ps.gfl(:,C.gfl.Sn),eps_v);
        Qref = ps.gfl(:,C.gfl.Qref) ./ max(ps.gfl(:,C.gfl.Sn),eps_v);
        id_ref = (2/3) .* Pref ./ vd_eff;
        iq_ref = -(2/3) .* Qref ./ vd_eff;
        Iref = sqrt(id_ref.^2 + iq_ref.^2);
        over_limit = Iref > ps.gfl(:,C.gfl.Imax);
        scale = ones(n_gfl,1);
        scale(over_limit) = ps.gfl(over_limit,C.gfl.Imax) ./ max(Iref(over_limit),eps_v);
        id_ref = id_ref .* scale;
        iq_ref = iq_ref .* scale;

        omega_hat = omega_0 + ps.gfl(:,C.gfl.Kp_pll).*vq + ps.gfl(:,C.gfl.Ki_pll).*xi_pll;
        f(ix.f.xi_pll_dot(gfl_on)) = vq(gfl_on);
        f(ix.f.rho_gfl_dot(gfl_on)) = omega_hat(gfl_on);
        f(ix.f.xi_id_dot(gfl_on)) = id_ref(gfl_on) - id(gfl_on);
        f(ix.f.xi_iq_dot(gfl_on)) = iq_ref(gfl_on) - iq(gfl_on);

        ud = ps.gfl(:,C.gfl.Kp_i).*(id_ref-id) + ps.gfl(:,C.gfl.Ki_i).*xi_id;
        uq = ps.gfl(:,C.gfl.Kp_i).*(iq_ref-iq) + ps.gfl(:,C.gfl.Ki_i).*xi_iq;
        Lf_eff = max(ps.gfl(:,C.gfl.Lf),eps_v);
        f(ix.f.id_gfl_dot(gfl_on)) = (-ps.gfl(gfl_on,C.gfl.Rf).*id(gfl_on) + ud(gfl_on)) ./ Lf_eff(gfl_on);
        f(ix.f.iq_gfl_dot(gfl_on)) = (-ps.gfl(gfl_on,C.gfl.Rf).*iq(gfl_on) + uq(gfl_on)) ./ Lf_eff(gfl_on);
    end
end

% output df_dx and df_dy if requested
if nargout>1
    % build df_dx
    if ~angle_ref
        dPg_ddelta = (Eaps.*mac_Vmags./Xdps).*cos(delta_m) + mac_Vmags.^2.*(1./Xqs-1./Xdps).*cos(2*delta_m);
        dFswing_ddelta_values = -dPg_ddelta./Ms;
        dFswing_domega_values = -Ds./Ms;
        dFdelta_dot_domega    = omega_0;
        dEap_ddelta_values    = -(Xds./Xdps-1).*mac_Vmags.*sin(delta_m)./Td0ps; 
    else
        dPg_ddelta = ((Eaps.*mac_Vmags./Xdps).*cos(delta_m) + mac_Vmags.^2.*(1./Xqs-1./Xdps).*cos(2*delta_m))*ones(1,ng).*((eye(ng)-weight./sum(weight)*ones(1,ng))');
        dFswing_ddelta_values = -dPg_ddelta./(Ms*ones(1,ng));
        dFswing_ddelta_values = dFswing_ddelta_values';
        dFswing_ddelta_values = dFswing_ddelta_values(:);
        
        dFswing_domega_values = -Ds*ones(1,ng).*(eye(ng)-weight./sum(weight)*ones(1,ng))'./(Ms*ones(1,ng));
        dFswing_domega_values = dFswing_domega_values';
        dFswing_domega_values = dFswing_domega_values(:);
        
        dFdelta_dot_domega    = omega_0.*(eye(ng)-weight./sum(weight)*ones(1,ng))';
        dFdelta_dot_domega    = dFdelta_dot_domega';
        dFdelta_dot_domega    = dFdelta_dot_domega(:);
        
        dEap_ddelta_values    = -(Xds./Xdps-1).*mac_Vmags.*sin(delta_m)./Td0ps*ones(1,ng).*((eye(ng)-weight./sum(weight)*ones(1,ng))');
        dEap_ddelta_values    = dEap_ddelta_values';
        dEap_ddelta_values    = dEap_ddelta_values(:);
    end

    dFswing_dPm_values    = 1./Ms;
    dFswing_dEa_values    = -sin(delta_m).*mac_Vmags./(Ms.*Xdps);
    dEap_dEap_values      = -Xds./(Td0ps.*Xdps);    
    dEap_dEfd_values      = 1./Td0ps;
    if gc
        dE1_dE1_values         = df_dx_exc(:,1);
        dEfd_dE1_values        = df_dx_exc(:,2);
        dEfd_dEfd_values       = df_dx_exc(:,3);
        dPm_domegas_values     = df_dx_gov(:,1);
        dPm_dPm_values         = df_dx_gov(:,2);
        dPm_dP3_values         = df_dx_gov(:,3);
        dP3_dP3_values         = df_dx_gov(:,4);
        dP3_domegas_values     = df_dx_gov(:,5);
    end

    if ~angle_ref
        omega_dot_loc = ix.f.omega_dot;
        delta_loc = ix.x.delta;
        delta_dot_loc = ix.f.delta_dot;
        omega_pu_loc = ix.x.omega_pu;
        Eap_dot_loc = ix.f.Eap_dot;
    else
        omega_dot_loc = ix.COI.omega_dot;
        delta_loc = ix.COI.delta;
        delta_dot_loc = ix.COI.delta_dot;
        omega_pu_loc = ix.COI.omega_pu;
        Eap_dot_loc = ix.COI.Eap_dot;
    end
    % assemble df_dx
    df_dx = sparse(ix.nx,ix.nx);
    % dFswing_ddelta
    df_dx = df_dx + sparse(omega_dot_loc,delta_loc,dFswing_ddelta_values,ix.nx,ix.nx);
    % dFswing_domega
    df_dx = df_dx + sparse(omega_dot_loc,omega_pu_loc,dFswing_domega_values,ix.nx,ix.nx);
    % dFswing_dPm
    df_dx = df_dx + sparse(ix.f.omega_dot,ix.x.Pm,dFswing_dPm_values,ix.nx,ix.nx);
    % dFswing_dEa
    df_dx = df_dx + sparse(ix.f.omega_dot,ix.x.Eap,dFswing_dEa_values,ix.nx,ix.nx);
    % dFdelta_dot_domega
    df_dx = df_dx + sparse(delta_dot_loc,omega_pu_loc,dFdelta_dot_domega,ix.nx,ix.nx);
	% dEap_dot_dEap
    df_dx = df_dx + sparse(ix.f.Eap_dot,ix.x.Eap,dEap_dEap_values,ix.nx,ix.nx);
    % dEap_dot_ddelta
    df_dx = df_dx + sparse(Eap_dot_loc,delta_loc,dEap_ddelta_values,ix.nx,ix.nx);
    % dEap_dot_dEfd
    df_dx = df_dx + sparse(ix.f.Eap_dot,ix.x.Efd,dEap_dEfd_values,ix.nx,ix.nx);
    if gc
        % dE1_dot_dE1
        df_dx = df_dx + sparse(ix.f.E1_dot,ix.x.E1,dE1_dE1_values,ix.nx,ix.nx);
        % dEfd_dot_dE1
        df_dx = df_dx + sparse(ix.f.Efd_dot,ix.x.E1,dEfd_dE1_values,ix.nx,ix.nx);
        % dEfd_dot_dEfd
        df_dx = df_dx + sparse(ix.f.Efd_dot,ix.x.Efd,dEfd_dEfd_values,ix.nx,ix.nx);
        % dPm_domega
        df_dx = df_dx + sparse(ix.f.Pm_dot,ix.x.omega_pu,dPm_domegas_values,ix.nx,ix.nx);
        % dPm_dot_dPm
        df_dx = df_dx + sparse(ix.f.Pm_dot,ix.x.Pm,dPm_dPm_values,ix.nx,ix.nx);
        % dPm_dP3
        df_dx = df_dx + sparse(ix.f.Pm_dot,ix.x.P3,dPm_dP3_values,ix.nx,ix.nx);
        % dP3_dot_dP3
        df_dx = df_dx + sparse(ix.f.P3_dot,ix.x.P3,dP3_dP3_values,ix.nx,ix.nx);
        % dP3_dot_domegas
        df_dx = df_dx + sparse(ix.f.P3_dot,ix.x.omega_pu,dP3_domegas_values,ix.nx,ix.nx);
    end
    % --- GFL df_dx dominant terms (new) ---
    if n_gfl > 0
        gfl_on = ps.gfl(:,C.gfl.status) > 0;
        if any(gfl_on)
            eps_v = 1e-4;
            rho = x(ix.x.rho_gfl);
            xi_pll = x(ix.x.xi_pll);
            id = x(ix.x.id_gfl);
            iq = x(ix.x.iq_gfl);
            gfl_bus_i = ps.bus_i(ps.gfl(:,C.gfl.bus));
            Vgfl = Vmags(gfl_bus_i);
            Tgfl = Thetas(gfl_bus_i);
            delta_theta = Tgfl - rho;
            vd = Vgfl .* cos(delta_theta);
            vq = Vgfl .* sin(delta_theta);
            vd_eff = max(vd,eps_v);
            Pref = ps.gfl(:,C.gfl.Pref) ./ max(ps.gfl(:,C.gfl.Sn),eps_v);
            Qref = ps.gfl(:,C.gfl.Qref) ./ max(ps.gfl(:,C.gfl.Sn),eps_v);
            id_ref = (2/3) .* Pref ./ vd_eff;
            iq_ref = -(2/3) .* Qref ./ vd_eff;
            Iref = sqrt(id_ref.^2 + iq_ref.^2);
            over_limit = Iref > ps.gfl(:,C.gfl.Imax);
            scale = ones(n_gfl,1);
            scale(over_limit) = ps.gfl(over_limit,C.gfl.Imax) ./ max(Iref(over_limit),eps_v);
            id_ref = id_ref .* scale;
            iq_ref = iq_ref .* scale;

            d_vq_d_rho = -Vgfl .* cos(delta_theta);
            d_vd_d_rho = Vgfl .* sin(delta_theta);
            d_idref_d_rho = -(2/3) .* Pref .* d_vd_d_rho ./ (vd_eff.^2);
            d_iqref_d_rho = +(2/3) .* Qref .* d_vd_d_rho ./ (vd_eff.^2);

            Lf_eff = max(ps.gfl(:,C.gfl.Lf),eps_v);
            d_did_dot_did = -(ps.gfl(:,C.gfl.Rf) + ps.gfl(:,C.gfl.Kp_i)) ./ Lf_eff;
            d_diq_dot_diq = -(ps.gfl(:,C.gfl.Rf) + ps.gfl(:,C.gfl.Kp_i)) ./ Lf_eff;
            d_did_dot_dxi_id = ps.gfl(:,C.gfl.Ki_i) ./ Lf_eff;
            d_diq_dot_dxi_iq = ps.gfl(:,C.gfl.Ki_i) ./ Lf_eff;
            d_did_dot_drho = ps.gfl(:,C.gfl.Kp_i) .* d_idref_d_rho ./ Lf_eff;
            d_diq_dot_drho = ps.gfl(:,C.gfl.Kp_i) .* d_iqref_d_rho ./ Lf_eff;

            df_dx = df_dx + sparse(ix.f.xi_pll_dot(gfl_on),ix.x.rho_gfl(gfl_on),d_vq_d_rho(gfl_on),ix.nx,ix.nx);
            df_dx = df_dx + sparse(ix.f.rho_gfl_dot(gfl_on),ix.x.rho_gfl(gfl_on),ps.gfl(gfl_on,C.gfl.Kp_pll).*d_vq_d_rho(gfl_on),ix.nx,ix.nx);
            df_dx = df_dx + sparse(ix.f.rho_gfl_dot(gfl_on),ix.x.xi_pll(gfl_on),ps.gfl(gfl_on,C.gfl.Ki_pll),ix.nx,ix.nx);
            df_dx = df_dx + sparse(ix.f.xi_id_dot(gfl_on),ix.x.id_gfl(gfl_on),-ones(sum(gfl_on),1),ix.nx,ix.nx);
            df_dx = df_dx + sparse(ix.f.xi_iq_dot(gfl_on),ix.x.iq_gfl(gfl_on),-ones(sum(gfl_on),1),ix.nx,ix.nx);
            df_dx = df_dx + sparse(ix.f.xi_id_dot(gfl_on),ix.x.rho_gfl(gfl_on),d_idref_d_rho(gfl_on),ix.nx,ix.nx);
            df_dx = df_dx + sparse(ix.f.xi_iq_dot(gfl_on),ix.x.rho_gfl(gfl_on),d_iqref_d_rho(gfl_on),ix.nx,ix.nx);
            df_dx = df_dx + sparse(ix.f.id_gfl_dot(gfl_on),ix.x.id_gfl(gfl_on),d_did_dot_did(gfl_on),ix.nx,ix.nx);
            df_dx = df_dx + sparse(ix.f.id_gfl_dot(gfl_on),ix.x.xi_id(gfl_on),d_did_dot_dxi_id(gfl_on),ix.nx,ix.nx);
            df_dx = df_dx + sparse(ix.f.id_gfl_dot(gfl_on),ix.x.rho_gfl(gfl_on),d_did_dot_drho(gfl_on),ix.nx,ix.nx);
            df_dx = df_dx + sparse(ix.f.iq_gfl_dot(gfl_on),ix.x.iq_gfl(gfl_on),d_diq_dot_diq(gfl_on),ix.nx,ix.nx);
            df_dx = df_dx + sparse(ix.f.iq_gfl_dot(gfl_on),ix.x.xi_iq(gfl_on),d_diq_dot_dxi_iq(gfl_on),ix.nx,ix.nx);
            df_dx = df_dx + sparse(ix.f.iq_gfl_dot(gfl_on),ix.x.rho_gfl(gfl_on),d_diq_dot_drho(gfl_on),ix.nx,ix.nx);
            % TODO: include derivatives through current limiting scaling when active.
        end
    end
end
if nargout>2
    % build df_dy (change in f wrt the algebraic variables)
    dPg_dVmag = Eaps.*sin(delta_m)./Xdps + mac_Vmags.*(1./Xqs-1./Xdps).*sin(2*delta_m);
    dFswing_dVmag_values        = -(dPg_dVmag)./Ms;
	dEap_dVmag_values           = (Xds./Xdps-1).*cos(delta_m)./Td0ps;
    if gc
        dE1_dVmag_values            = df_dy_exc(:,1);
        dEfd_dVmag_values           = df_dy_exc(:,2);
    end
    
	% assemble df_dy
    cols = ix.y.Vmag(mac_bus_i);
    df_dy = sparse(ix.f.omega_dot,cols,dFswing_dVmag_values,ix.nx,ix.ny);
    
    if ~angle_ref
        dFswing_dtheta_values = dFswing_ddelta_values;
        cols = ix.y.delta_sys;
        % domega_dot_ddelta_sys
        df_dy = df_dy + sparse(ix.f.omega_dot,cols,dFswing_ddelta_values,ix.nx,ix.ny);
        %dEap_dot_ddelta_sys
        cols = ix.y.delta_sys;
        df_dy = df_dy + sparse(ix.f.Eap_dot,cols,dEap_ddelta_values,ix.nx,ix.ny);
    else
        dFswing_dtheta_values = -((Eaps.*mac_Vmags./Xdps).*cos(delta_m) + mac_Vmags.^2.*(1./Xqs-1./Xdps).*cos(2*delta_m))./Ms;
    end
    cols = ix.y.theta(mac_bus_i);
    df_dy = df_dy + sparse(ix.f.omega_dot,cols,-dFswing_dtheta_values,ix.nx,ix.ny);
	
	% dEap_dot_dVmag
    cols  = ix.y.Vmag(mac_bus_i);
    df_dy = df_dy + sparse(ix.f.Eap_dot,cols,dEap_dVmag_values,ix.nx,ix.ny);
    % dEap_dot_dtheta
    cols = ix.y.theta(mac_bus_i);
    dEap_dtheta_values = -(Xds./Xdps-1).*mac_Vmags.*sin(delta_m)./Td0ps;
    df_dy = df_dy + sparse(ix.f.Eap_dot,cols,-dEap_dtheta_values,ix.nx,ix.ny);

    if gc
        % dE1_dot_dVmag
        cols  = ix.y.Vmag(mac_bus_i);
        df_dy = df_dy + sparse(ix.f.E1_dot,cols,dE1_dVmag_values,ix.nx,ix.ny);
        % dEfd_dot_dVmag
        cols  = ix.y.Vmag(mac_bus_i);
        df_dy = df_dy + sparse(ix.f.Efd_dot,cols,dEfd_dVmag_values,ix.nx,ix.ny);
    end
    % --- GFL df_dy dominant terms (new) ---
    if n_gfl > 0
        gfl_on = ps.gfl(:,C.gfl.status) > 0;
        if any(gfl_on)
            eps_v = 1e-4;
            rho = x(ix.x.rho_gfl);
            xi_pll = x(ix.x.xi_pll); %#ok<NASGU>
            id = x(ix.x.id_gfl); %#ok<NASGU>
            iq = x(ix.x.iq_gfl); %#ok<NASGU>
            gfl_bus_i = ps.bus_i(ps.gfl(:,C.gfl.bus));
            Vgfl = Vmags(gfl_bus_i);
            Tgfl = Thetas(gfl_bus_i);
            delta_theta = Tgfl - rho;
            vd = Vgfl .* cos(delta_theta);
            vq = Vgfl .* sin(delta_theta);
            vd_eff = max(vd,eps_v);
            Pref = ps.gfl(:,C.gfl.Pref) ./ max(ps.gfl(:,C.gfl.Sn),eps_v);
            Qref = ps.gfl(:,C.gfl.Qref) ./ max(ps.gfl(:,C.gfl.Sn),eps_v);
            d_vq_dV = sin(delta_theta);
            d_vq_dTheta = Vgfl .* cos(delta_theta);
            d_vd_dV = cos(delta_theta);
            d_vd_dTheta = -Vgfl .* sin(delta_theta);
            d_idref_dV = -(2/3) .* Pref .* d_vd_dV ./ (vd_eff.^2);
            d_iqref_dV = +(2/3) .* Qref .* d_vd_dV ./ (vd_eff.^2);
            d_idref_dTheta = -(2/3) .* Pref .* d_vd_dTheta ./ (vd_eff.^2);
            d_iqref_dTheta = +(2/3) .* Qref .* d_vd_dTheta ./ (vd_eff.^2);
            Lf_eff = max(ps.gfl(:,C.gfl.Lf),eps_v);
            d_did_dot_dV = ps.gfl(:,C.gfl.Kp_i) .* d_idref_dV ./ Lf_eff;
            d_diq_dot_dV = ps.gfl(:,C.gfl.Kp_i) .* d_iqref_dV ./ Lf_eff;
            d_did_dot_dTheta = ps.gfl(:,C.gfl.Kp_i) .* d_idref_dTheta ./ Lf_eff;
            d_diq_dot_dTheta = ps.gfl(:,C.gfl.Kp_i) .* d_iqref_dTheta ./ Lf_eff;

            colsV = ix.y.Vmag(gfl_bus_i(gfl_on));
            colsT = ix.y.theta(gfl_bus_i(gfl_on));
            df_dy = df_dy + sparse(ix.f.xi_pll_dot(gfl_on),colsV,d_vq_dV(gfl_on),ix.nx,ix.ny);
            df_dy = df_dy + sparse(ix.f.xi_pll_dot(gfl_on),colsT,d_vq_dTheta(gfl_on),ix.nx,ix.ny);
            df_dy = df_dy + sparse(ix.f.rho_gfl_dot(gfl_on),colsV,ps.gfl(gfl_on,C.gfl.Kp_pll).*d_vq_dV(gfl_on),ix.nx,ix.ny);
            df_dy = df_dy + sparse(ix.f.rho_gfl_dot(gfl_on),colsT,ps.gfl(gfl_on,C.gfl.Kp_pll).*d_vq_dTheta(gfl_on),ix.nx,ix.ny);
            df_dy = df_dy + sparse(ix.f.xi_id_dot(gfl_on),colsV,d_idref_dV(gfl_on),ix.nx,ix.ny);
            df_dy = df_dy + sparse(ix.f.xi_id_dot(gfl_on),colsT,d_idref_dTheta(gfl_on),ix.nx,ix.ny);
            df_dy = df_dy + sparse(ix.f.xi_iq_dot(gfl_on),colsV,d_iqref_dV(gfl_on),ix.nx,ix.ny);
            df_dy = df_dy + sparse(ix.f.xi_iq_dot(gfl_on),colsT,d_iqref_dTheta(gfl_on),ix.nx,ix.ny);
            df_dy = df_dy + sparse(ix.f.id_gfl_dot(gfl_on),colsV,d_did_dot_dV(gfl_on),ix.nx,ix.ny);
            df_dy = df_dy + sparse(ix.f.id_gfl_dot(gfl_on),colsT,d_did_dot_dTheta(gfl_on),ix.nx,ix.ny);
            df_dy = df_dy + sparse(ix.f.iq_gfl_dot(gfl_on),colsV,d_diq_dot_dV(gfl_on),ix.nx,ix.ny);
            df_dy = df_dy + sparse(ix.f.iq_gfl_dot(gfl_on),colsT,d_diq_dot_dTheta(gfl_on),ix.nx,ix.ny);
            % TODO: include df_dy terms from current limit scaling.
        end
    end
end
