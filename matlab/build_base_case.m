function [ps, opt] = build_base_case()
    C = psconstants;

    ps = updateps(case39_ps);
    ps = unify_generators(ps);

    % 基础设置
    ps.branch(:,C.br.tap)   = 1;
    ps.shunt(:,C.sh.factor) = 1;
    ps.shunt(:,C.sh.status) = 1;

    opt = psoptions;
    opt.sim.integration_scheme = 1;
    opt.sim.dt_default = 0.1;
    opt.verbose = false;
    opt.sim.gen_control = 1;

    % 保护延时
    opt.sim.uvls_tdelay_ini = 0.5;
    opt.sim.ufls_tdelay_ini = 0.5;
    opt.sim.dist_tdelay_ini = 0.5;
    opt.sim.temp_tdelay_ini = 0;

    % 初始化潮流和动态对象
    ps = newpf(ps,opt);
    [ps.Ybus,ps.Yf,ps.Yt] = getYbus(ps,false);
    ps = update_load_freq_source(ps);
    [ps.mac,ps.exc,ps.gov] = get_mac_state(ps,'salient');
    ps.relay = get_relays(ps,'all',opt);
end