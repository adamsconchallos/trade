function [tFH_N, tHF_N, info] = nash_tariffs_best_response(par, t0, tmax, tol, maxit, wH_guess)
    if nargin < 2, t0   = 0.05; end
    if nargin < 3, tmax = 2.0;  end
    if nargin < 4, tol  = 1e-6; end
    if nargin < 5, maxit= 200;  end
    if nargin < 6, wH_guess = 1.0; end

    tFH = t0; tHF = t0;
    for it = 1:maxit
        % Home best response to current Foreign tariff
        [tFH_new, ~] = opt_tariff_H(par, tHF, tmax, wH_guess);

        % Foreign best response: swap country labels logically by symmetry.
        % Implement by *reusing* the same routines after swapping H <-> F.
        parF = par;  % build a relabeled problem
        parF.aH = par.aF; parF.aF = par.aH;
        parF.LH = par.LF; parF.LF = par.LH;
        % Note: when we solve Foreign's problem, "Home" in the helper is actually Foreign.
        [tHF_new, ~] = opt_tariff_H(parF, tFH_new, tmax, wH_guess);

        % Convergence check
        if max(abs([tFH_new - tFH, tHF_new - tHF])) < tol
            tFH = tFH_new; tHF = tHF_new;
            info = struct('iter', it, 'converged', true);
            tFH_N = tFH; tHF_N = tHF; return;
        end
        tFH = tFH_new; tHF = tHF_new;
    end
    info = struct('iter', maxit, 'converged', false);
    tFH_N = tFH; tHF_N = tHF;
end
