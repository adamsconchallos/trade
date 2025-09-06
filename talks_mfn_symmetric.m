function [t_coop, WH, WF, sol] = talks_mfn_symmetric(par, tmax, wH_guess)
    if nargin < 2, tmax = 2.0; end
    if nargin < 3, wH_guess = 1.0; end

    % Threat point: Nash (trade war)
    [tFH_N, tHF_N] = nash_tariffs_best_response(par);
    [WH_N, WF_N, ~] = welfare_given_tariffs(par, tFH_N, tHF_N, wH_guess);

    % Nash–bargaining objective with MFN+reciprocity: tFH = tHF = t
    obj = @(t) -barg_obj_sym(t, par, WH_N, WF_N, wH_guess);

    % Solve on [0, tmax]
    [t_coop, ~] = fminbnd(obj, 0, tmax);

    % Cooperative outcome
    [WH, WF, sol] = welfare_given_tariffs(par, t_coop, t_coop, wH_guess);
end

function val = barg_obj_sym(t, par, WH_N, WF_N, wH_guess)
    [WH, WF] = welfare_given_tariffs(par, t, t, wH_guess);
    gH = WH - WH_N;  gF = WF - WF_N;
    if gH <= 0 || gF <= 0
        val = 1e6;                          % outside bargaining set
    else
        val = -(gH * gF);                   % maximize product of gains
    end
end
