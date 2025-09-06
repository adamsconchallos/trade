function [WH, WF, sol] = welfare_given_tariffs(par, tFH, tHF, wH_guess)
    par.tFH = tFH; par.tHF = tHF;
    if nargin < 4, wH_guess = 1.0; end
    resfun = @(w) ge_armington_two_country(w, par);
    wH = fzero(resfun, wH_guess);
    [~, out] = ge_armington_two_country(wH, par);
    WH = out.W.H; WF = out.W.F;
    sol = struct('wH', wH, 'out', out);
end
