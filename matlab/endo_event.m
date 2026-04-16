function [value,isterminal,direction,Temperature] = endo_event(t,xy,ix,ps,dt,opt)
% usage: [value,isterminal,direction] = endo_event(t,xy,ix,ps)
% defines an exogenous event (relay) that will stop the DAE integration

% constants and settings
global dist2threshold state_a temp gfl_rocof_state
C = psconstants;
oc_setting1             = ps.relay(ix.re.oc,C.re.setting1);
Vmag_threshold          = ps.relay(ix.re.uvls,C.re.threshold);
omega_pu_threshold      = ps.relay(ix.re.ufls,C.re.threshold);
dist_zone1_threshold    = ps.relay(ix.re.dist,C.re.threshold);
temp_threshold          = ps.relay(ix.re.temp,C.re.threshold);
temp_K                  = ps.relay(ix.re.temp,C.re.temp_K);
temp_R                  = ps.relay(ix.re.temp,C.re.temp_R);
SMALL_EPS = 1e-12;
omega_0 = 2*pi*ps.frequency;

x               = xy(1:ix.nx);
y               = xy(ix.nx+1:end);
Vmags           = y(ix.y.Vmag);
Thetas          = y(ix.y.theta);
V               = Vmags .* exp(1i*Thetas);
If              = ps.Yf * V;
It              = ps.Yt * V;
Imag_f          = abs(If);
Imag_t          = abs(It);
Imag            = max(Imag_f, Imag_t);

temp(ps.relay(ix.re.temp,C.re.id),1)  = x(ix.x.temp);
temp(ps.relay(ps.branch(:,C.br.status)==0,C.re.id),1) = 0;

if ~isempty(dt)
    state_a(ps.relay(ix.re.oc,C.re.id)) = max(state_a(ps.relay(ix.re.oc,C.re.id))+(Imag-oc_setting1)*dt+SMALL_EPS, 0);
    temp(ps.relay(ix.re.temp,C.re.id))  = temp(ps.relay(ix.re.temp,C.re.id)) + (temp_R.*Imag_f.^2 -temp_K.*temp(ps.relay(ix.re.temp,C.re.id)))*dt;
end
Temperature = temp(ps.relay(ix.re.temp,C.re.id));

dist2threshold(ps.relay(ix.re.oc,C.re.id)) = ps.relay(ix.re.oc,C.re.threshold) - state_a(ps.relay(ix.re.oc,C.re.id));
dist2threshold(dist2threshold<0)=0;

F               = ps.bus_i(ps.branch(:,C.br.from));
y_apparent      = Imag_f./Vmags(F);
nload           = size(ps.shunt,1);
load_freq       = zeros(nload,1);
near_gen        = ps.shunt(:,C.sh.near_gen);
for i = 1:nload
    near_gen_id = ps.gen_i(near_gen(i));
    if near_gen_id ~= 0
        load_freq_source = ix.x.omega_pu(near_gen_id);
        load_freq(i)  = x(load_freq_source);
    end
end

sh_bus_ix = ps.bus_i(ps.shunt(:,1));
Vmag_sh   = y(ix.y.Vmag(sh_bus_ix));

value = [temp_threshold - Temperature;
 oc_setting1          - Imag;
 Vmag_sh              - Vmag_threshold;
 load_freq            - omega_pu_threshold;
 dist_zone1_threshold - y_apparent];

if isfield(ps,'gfl') && ~isempty(ps.gfl) && isfield(ix.re,'gfl_rocof') && ~isempty(ix.re.gfl_rocof) ...
        && size(ps.relay,1) >= ix.re.gfl_rocof(end)
    n_gfl = size(ps.gfl,1);
    gfl_on = ps.gfl(:,C.gfl.status) > 0;
    if n_gfl > 0
        rho = x(ix.x.rho_gfl);
        xi_pll = x(ix.x.xi_pll);
        id = x(ix.x.id_gfl);
        iq = x(ix.x.iq_gfl);
        gfl_bus_i = ps.bus_i(ps.gfl(:,C.gfl.bus));
        Vgfl = Vmags(gfl_bus_i);
        Tgfl = Thetas(gfl_bus_i);
        vq_gfl = Vgfl .* sin(Tgfl-rho);
        omega_hat = omega_0 + ps.gfl(:,C.gfl.Kp_pll).*vq_gfl + ps.gfl(:,C.gfl.Ki_pll).*xi_pll;
        pll_dev = abs(omega_hat-omega_0);
        Igfl = sqrt(id.^2 + iq.^2);
        val_uv = Vgfl - ps.gfl(:,C.gfl.Vtrip);
        val_oc = ps.gfl(:,C.gfl.Imax) - Igfl;
        val_pll = ps.gfl(:,C.gfl.wmax_dev) - pll_dev;

        rocof_abs = zeros(n_gfl,1);
        if isempty(gfl_rocof_state) || ~isfield(gfl_rocof_state,'omega_hat') || numel(gfl_rocof_state.omega_hat) ~= n_gfl
            gfl_rocof_state.omega_hat = omega_hat;
            gfl_rocof_state.t = t;
        else
            dt_rocof = dt;
            if isempty(dt_rocof)
                dt_rocof = t - gfl_rocof_state.t;
            end
            if ~isempty(dt_rocof) && dt_rocof > 0
                rocof_abs = abs((omega_hat - gfl_rocof_state.omega_hat) ./ dt_rocof);
            end
            if ~isempty(dt) && dt > 0
                gfl_rocof_state.omega_hat = omega_hat;
                gfl_rocof_state.t = t;
            end
        end

        rocof_ids = ps.relay(ix.re.gfl_rocof,C.re.id);
        max_rocof_id = max(rocof_ids);
        if numel(state_a) < max_rocof_id
            state_a(max_rocof_id,1) = 0;
        end
        if numel(dist2threshold) < max_rocof_id
            dist2threshold(max_rocof_id,1) = inf;
        end
        rocof_latched = max(state_a(rocof_ids),rocof_abs);
        if ~isempty(dt) && dt > 0
            state_a(rocof_ids) = rocof_latched;
            dist2threshold(rocof_ids) = ps.relay(ix.re.gfl_rocof,C.re.threshold) - rocof_latched;
        end
        val_rocof = ps.relay(ix.re.gfl_rocof,C.re.threshold) - rocof_latched;

        val_uv(~gfl_on) = 10;
        val_oc(~gfl_on) = 10;
        val_pll(~gfl_on) = 10;
        val_rocof(~gfl_on) = 10;
        value = [value; val_uv; val_oc; val_pll; val_rocof];
    end
end

is_tripped = ps.relay(:,C.re.tripped)==1;
value(is_tripped) = 10;

n_relays = size(ps.relay,1);
isterminal  = ones(n_relays,1);
direction   = -ones(n_relays,1);
