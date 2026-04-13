function cont = generate_random_contingency(ps, mode)
    C = psconstants;

    in_service = find(ps.branch(:,C.br.status) == 1);
    branch_ids = ps.branch(in_service, C.br.id);

    if strcmpi(mode,'N-1')
        k = 1;
    else
        k = 2;
    end

    pick = randperm(length(branch_ids), k);

    cont.type = mode;
    cont.branch_ids = branch_ids(pick);
    cont.t_fault = 3.0;
end