%----------------------------------------------------
% Constraints for cooperative (non-MFN) tariff problem
%   x = [ G_hat (N) ; X_hat (N) ; W_hat (N) ; C_transf ((N-1)*N*S) ]
% where C_transf are transformed off-diagonal import tariffs.
%----------------------------------------------------
function [c, ceq] = myconcoop(x)

    % Globals:
    %  - N, S      from mycalculations
    %  - ARG       so that LAMBDA = ARG + 1 (political weights)
    %  - SCALING   tariff scaling constant
    %  - TAU       scaling tensor for transformed tariffs, size (N-1)×N×S
    % If your solver still defines TAUs (plural), change TAU→TAUs below.
    global N S
    global ARG SCALING TAU

    % ---- Unpack decision vector ----
    LAMBDA = ARG + 1;                      % N×S
    Ghat   = x(1:N);                       % government welfare (levels)
    Xhat   = x(N+1:2*N);                   % expenditure (levels)
    What   = x(2*N+1:3*N);                 % wages (levels)
    Ctrans = x(3*N+1:end);                 % transformed tariffs, vector

    % ---- Back-transform tariffs to levels on the off-diagonal stack ----
    % Ctrans has length (N-1)*N*S; reshape to (N-1)×N×S to align with TAU
    Cstack = reshape(Ctrans, [N-1, N, S]);             % (N-1)×N×S
    % If your globals use TAUs (plural), replace TAU with TAUs here:
    TEMP   = (Cstack .* TAU - 1) ./ SCALING;           % (N-1)×N×S

    % ---- Inflate to bilateral cube with zero diagonals: N×N×S ----
    TARIFFCs = zeros(N, N, S);
    for j = 1:N
        TARIFFCs(:, j, :) = [TEMP(1:j-1, j, :); zeros(1,1,S); TEMP(j:end, j, :)];
    end

    % ---- GE conditions given (Xhat, What, TARIFFCs) ----
    % N.B.: myconditions returns Walras + numeraire; we drop one equation
    [EQ, NUM, ~, PaggHat, TRADECs] = myconditions(Xhat, What, TARIFFCs, zeros(N,1));
    EQ(1) = [];                               % drop one equilibrium condition (Walras)

    % ---- Government welfare consistency ----
    Ghat_cons = mycongov(LAMBDA, TRADECs, TARIFFCs, PaggHat);

    % ---- Equalities: market-clearing/numeraire + gov’t welfare + normalization ----
    ceq1 = [EQ; NUM];                         % (2N-1)×1
    ceq2 = Ghat - Ghat_cons;                  % N×1

    % Normalization/identification on government welfare (original code):
    % [1, -I] * Ghat enforces G1 - G2..N = 0 → pins levels relative to country 1
    ceq3 = [ones(N-1,1) -eye(N-1)] * Ghat;    % (N-1)×1

    % ---- Outputs ----
    c   = [];                                 % no inequalities
    ceq = [ceq1; ceq2; ceq3];

    % ---- Optional size checks (uncomment if debugging) ----
    % assert(all(size(TAU) == [N-1, N, S]), 'TAU must be (N-1)×N×S');
    % assert(numel(Ctrans) == (N-1)*N*S,    'C_transf length mismatch');

end
