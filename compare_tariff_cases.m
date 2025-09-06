function T = compare_tariff_cases(par, tbox)
% Returns a table with tariffs, terms of trade, welfare, and %Δ vs free trade.
% tbox is either a scalar upper bound or a 2x2 bounds matrix for talks_bargaining_twovar.

    if nargin < 2, tbox = 1.0; end

    % ---------- 1) Free trade ----------
    par.tFH = 0; par.tHF = 0;
    wH = fzero(@(w) ge_armington_two_country(w,par), 1.0);
    [~, outFT] = ge_armington_two_country(wH, par);
    ToT_FT = outFT.p.H/outFT.p.F;
    WH_FT  = outFT.W.H;  WF_FT = outFT.W.F;

    % ---------- 2) Unilateral (Home best response to tHF=0) ----------
    [tFH_star, WH_uni, solU] = opt_tariff_H(par, 0.0, 1.0, 1.0);
    outUNI = solU.out;  WF_uni = outUNI.W.F;  ToT_UNI = outUNI.p.H/outUNI.p.F;

    % ---------- 3) Nash (trade war) ----------
    [tFH_N, tHF_N, info] = nash_tariffs_best_response(par);
    [WH_N, WF_N, solN]  = welfare_given_tariffs(par, tFH_N, tHF_N, 1.0);
    outN = solN.out;  ToT_N = outN.p.H/outN.p.F;

    % ---------- 4) Talks (two-variable bargaining) ----------
    [tFH_C, tHF_C, WH_C, WF_C, solC] = talks_bargaining_twovar(par, tbox, 1.0);
    ToT_C = solC.out.p.H/solC.out.p.F;

    % ---------- 5) Talks (MFN symmetric) ----------
    [t_sym, WH_sym, WF_sym, solSym] = talks_mfn_symmetric(par, 1.0, 1.0);
    ToT_sym = solSym.out.p.H/solSym.out.p.F;

    % ---------- Build table ----------
    rows = {
      'Free trade',        0,        0,        ToT_FT,  WH_FT,   WF_FT;
      'Home unilateral',   tFH_star, 0,        ToT_UNI, WH_uni,  WF_uni;
      'Nash',              tFH_N,    tHF_N,    ToT_N,   WH_N,    WF_N;
      'Talks (2-var)',     tFH_C,    tHF_C,    ToT_C,   WH_C,    WF_C;
      'Talks (MFN sym)',   t_sym,    t_sym,    ToT_sym, WH_sym,  WF_sym
    };
    T = cell2table(rows, 'VariableNames', ...
        {'Case','tFH','tHF','ToT_border','W_H','W_F'});

    % ---------- Add % changes vs Free Trade ----------
    T.dWH_pct_vs_FT = 100*(T.W_H/WH_FT - 1);
    T.dWF_pct_vs_FT = 100*(T.W_F/WF_FT - 1);
end
