function [t_star, WH_star, WF_star, sol] = opt_tariff_symmetric(par, tmax, wH_guess)
    if nargin < 2, tmax = 2.0; end
    if nargin < 3, wH_guess = 1.0; end
    obj = @(t) -welfare_given_tariffs(par, t, t, wH_guess); % Home’s welfare under symmetric t
    [t_star, negWH] = fminbnd(@(t) obj(t), 0, tmax);
    [WH_star, WF_star, sol] = welfare_given_tariffs(par, t_star, t_star, wH_guess);
end
