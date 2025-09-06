function sol = solve_GE_NS(par, w0)
    if nargin<2, w0 = ones(par.N-1,1); end
    resfun = @(wf) ge_armington_NS(wf, par);
    opts = optimoptions('fsolve','Display','iter','FunctionTolerance',1e-11,'StepTolerance',1e-11);
    [w_free,~,exitflag] = fsolve(@(wf) resfun(wf), w0, opts);
    [~, out] = ge_armington_NS(w_free, par);
    sol = struct('exitflag',exitflag,'out',out);
end
