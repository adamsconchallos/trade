function [tFH_star, WH_star, sol] = opt_tariff_H(par, tHF_fixed, tmax, wH_guess)
    if nargin < 3, tmax = 2.0; end         % allow up to 200%
    if nargin < 4, wH_guess = 1.0; end
    obj = @(tFH) -welfare_given_tariffs(par, tFH, tHF_fixed, wH_guess);
    [tFH_star, negWH] = fminbnd(@(t) obj(t), 0, tmax);
    [WH_star, ~, sol] = welfare_given_tariffs(par, tFH_star, tHF_fixed, wH_guess);
end
