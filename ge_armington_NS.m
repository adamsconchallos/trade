function [res, out] = ge_armington_NS(w_free, par)
% Multi-country (N), multi-sector (S) Armington GE.
% Numeraire: w(N)=1. Input vector w_free = w(1:N-1).

    % Unpack
    N = par.N; S = par.S;
    a = par.a;            % N×S
    L = par.L(:);         % N×1
    sigma = par.sigma(:); % S×1
    mu = par.mu;          % N×S, rows sum to 1
    kappa = par.kappa;    % N×N×S
    tau = par.tau;        % N×N×S  (ad valorem tariffs)
    
    % Wages incl. numeraire
    w = ones(N,1); w(1:N-1) = w_free;
    
    % Producer prices p_{i s} = w_i * a_{i s}
    p = w .* a;           % N×S
    
    % Consumer prices p^c_{i j s} = (1+tau_{i j s}) * kappa_{i j s} * p_{i s}
    % Build as 3D arrays over (i,j,s)
    pc = zeros(N,N,S);
    for s = 1:S
        ps = p(:,s);                            % N×1
        pc(:,:,s) = (1 + tau(:,:,s)) .* kappa(:,:,s) .* (ps * ones(1,N)); % (i,j)
    end
    
    % CES denominators and shares π_{i j s}
    den = zeros(N,N,S);  pi = zeros(N,N,S);
    for s = 1:S
        tpow = pc(:,:,s).^(1 - sigma(s));  % N×N
        den(:,:,s) = ones(N,1) * sum(tpow,1); % 1×N replicated down
        pi(:,:,s)  = tpow ./ den(:,:,s);      % N×N
    end
    
    % National incomes with tariff rebates (closure)
    % X_j = w_j L_j / (1 - sum_s sum_i [ t/(1+t) * μ_{j s} * π_{i j s} ])
    X = zeros(N,1);
    for j = 1:N
        rebate_share = 0.0;
        for s = 1:S
            rebate_share = rebate_share + mu(j,s) * sum( (tau(:,j,s)./(1+tau(:,j,s))) .* pi(:,j,s) );
        end
        X(j) = (w(j) * L(j)) / (1 - rebate_share);
    end
    
    % Sectoral expenditures and quantities
    % q_{i j s} = (π_{i j s} * μ_{j s} * X_j) / [(1+tau)*kappa*p_{i s}]
    q = zeros(N,N,S);  y = zeros(N,S);
    for s = 1:S
        denom_ijs = pc(:,:,s);         % (1+tau)*kappa*p
        for j = 1:N
            q(:,j,s) = (pi(:,j,s) * (mu(j,s)*X(j))) ./ denom_ijs(:,j);
        end
        y(:,s) = sum(q(:,:,s), 2);     % output by origin in sector s
    end
    
    % Labor demands and residuals
    LD = sum(a .* y, 2);   % N×1
    res_full = LD - L;
    res = res_full(1:N-1); % drop numeraire country residual
    
    % Price indices and welfare
    Psec = zeros(N,S);
    for s = 1:S
        Psec(:,s) = (sum(pc(:,:,s).^(1 - sigma(s)), 1).').^(1/(1 - sigma(s))); % N×1
    end
    P = exp( (log(Psec) .* mu) * ones(S,1) );   % P_j = ∏_s P_{js}^{μ_js}
    W = X ./ P;

    % Pack diagnostics
    out = struct('w',w,'p',p,'pc',pc,'pi',pi,'X',X,'y',y,'LD',LD,'Psec',Psec,'P',P,'W',W);
end
