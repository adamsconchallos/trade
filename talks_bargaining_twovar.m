function [tFH_C, tHF_C, WH, WF, sol] = talks_bargaining_twovar(par, box, wH_guess)
    if nargin<2, box=[0 2; 0 2]; end
    if nargin<3, wH_guess=1.0; end

    % Threat point = Nash
    [tFH_N, tHF_N] = nash_tariffs_best_response(par);
    [WH_N, WF_N]   = welfare_given_tariffs(par, tFH_N, tHF_N, wH_guess);

    % Product-of-gains objective
    function f = obj(x)
        tFH=max(x(1),0); tHF=max(x(2),0);
        [WHc,WFc] = welfare_given_tariffs(par, tFH, tHF, wH_guess);
        gH=max(WHc-WH_N,0); gF=max(WFc-WF_N,0);
        if gH==0 || gF==0, f=1e6; else, f=-(gH*gF); end
    end

    % Use pattern search or fmincon; here: fmincon with simple bounds
    x0=[tFH_N, tHF_N];
    lb=box(:,1)'; ub=box(:,2)';
    opts=optimoptions('fmincon','Display','iter','OptimalityTolerance',1e-8);
    [x, ~] = fmincon(@obj, x0, [],[],[],[], lb, ub, [], opts);

    tFH_C = x(1); tHF_C = x(2);
    [WH, WF, sol] = welfare_given_tariffs(par, tFH_C, tHF_C, wH_guess);
end
