function [res, out] = ge_armington_two_country(wH, par)
    % Numeraire and parameters
    wF  = 1.0;
    aH  = par.aH;    aF  = par.aF;
    LH  = par.LH;    LF  = par.LF;
    s   = par.sigma;
    tHF = par.tHF;   tFH = par.tFH;

    % Producer prices
    pH = wH * aH;
    pF = aF;  % since wF = 1

    % Consumer prices
    pHH = pH;              pFH = (1 + tFH) * pF;
    pHF = (1 + tHF) * pH;  pFF = pF;

    % CES shares (Home)
    denH  = pHH^(1-s) + pFH^(1-s);
    pi_HH = pHH^(1-s) / denH;
    pi_FH = pFH^(1-s) / denH;
    % CES shares (Foreign)
    denF  = pHF^(1-s) + pFF^(1-s);
    pi_HF = pHF^(1-s) / denF;
    pi_FF = pFF^(1-s) / denF;

    % Incomes with tariff rebates: t/(1+t) * (import share) * X_j
    XH = (wH * LH) / (1 - (tFH/(1+tFH)) * pi_FH);
    XF = (wF * LF) / (1 - (tHF/(1+tHF)) * pi_HF);

    % Quantities and output (producer receives border price)
    q_HH = (pi_HH * XH) / pHH;
    q_HF = (pi_HF * XF) / ((1 + tHF) * pH);
    yH   = q_HH + q_HF;

    LDH  = aH * yH;
    res  = LDH - LH;

    % Foreign block (diagnostics)
    q_FF = (pi_FF * XF) / pFF;
    q_FH = (pi_FH * XH) / ((1 + tFH) * pF);
    yF   = q_FF + q_FH;  LDF = aF * yF;

    % CES price indices and welfare
    PH = (pHH^(1-s) + pFH^(1-s))^(1/(1-s));
    PF = (pHF^(1-s) + pFF^(1-s))^(1/(1-s));

    out.p  = struct('H',pH,'F',pF,'HH',pHH,'FH',pFH,'HF',pHF,'FF',pFF);
    out.X  = struct('H',XH,'F',XF);
    out.q  = struct('HH',q_HH,'HF',q_HF,'FF',q_FF,'FH',q_FH);
    out.L  = struct('useH',LDH,'useF',LDF);
    out.P  = struct('H',PH,'F',PF);
    out.W  = struct('H',XH/PH,'F',XF/PF);
    out.shares = struct('pi_HH',pi_HH,'pi_FH',pi_FH,'pi_HF',pi_HF,'pi_FF',pi_FF);
end
