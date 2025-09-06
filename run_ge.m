% Baseline parameters
par.sigma = 4;
par.aH = 1.0; par.aF = 1.2;
par.LH = 1.0; par.LF = 1.0;

% ----- 1) Free trade -----
par.tFH = 0; par.tHF = 0;
wH = fzero(@(w) ge_armington_two_country(w,par), 1.0);
[~, outFT] = ge_armington_two_country(wH, par);

% ----- 2) Home unilateral optimum (Foreign tariff = 0) -----
[tFH_star, WH_star, solU] = opt_tariff_H(par, 0.0, 1.0);
outUNI = solU.out;
WF_star = outUNI.W.F;  % to see what happens to Foreign

% ----- 3) Nash equilibrium -----
[tFH_N, tHF_N, info] = nash_tariffs_best_response(par);
[WH_N, WF_N, solN] = welfare_given_tariffs(par, tFH_N, tHF_N, 1.0);
outN = solN.out;

% --- Talks (two-variable bargaining) ---
[tFH_C, tHF_C, WH_C, WF_C, solC] = talks_bargaining_twovar(par, [0 1; 0 1], 1.0);

% --- Talks (MFN symmetric) ---
[t_sym, WH_sym, WF_sym, solSym] = talks_mfn_symmetric(par, 1.0, 1.0);

% ----- Collect results (now including symmetric talks) -----
results = {
    'Free trade',        0,        0,        outFT.p.H/outFT.p.F,      outFT.W.H, outFT.W.F;
    'Home unilateral',   tFH_star, 0,        outUNI.p.H/outUNI.p.F,     outUNI.W.H, outUNI.W.F;
    'Nash',              tFH_N,    tHF_N,    outN.p.H/outN.p.F,         WH_N,       WF_N;
    'Talks (2-var)',     tFH_C,    tHF_C,    solC.out.p.H/solC.out.p.F, WH_C,       WF_C;
    'Talks (MFN sym)',   t_sym,    t_sym,    solSym.out.p.H/solSym.out.p.F, WH_sym, WF_sym
};

T = cell2table(results, 'VariableNames', ...
    {'Case','tFH','tHF','ToT_border','W_H','W_F'});
disp(T)
