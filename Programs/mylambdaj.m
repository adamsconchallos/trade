function [lambda_j, mfnTariff_j] = mylambdaj(j, targetTariff_j)
%MYLAMBDaj Compute political-economy parameter and MFN optimal tariff for country j
%
%   [lambda_j, mfnTariff_j] = mylambdaj(j, targetTariff_j)
%
%   Inputs:
%     j                 Index of the home country (1..N)
%     targetTariff_j    [S×1] vector of target ad-valorem tariffs for country j
%
%   Outputs:
%     lambda_j          [S×1] vector of political-economy parameters
%     mfnTariff_j       [S×1] vector of MFN-optimal tariffs for country j

  %% 1) Load data and dimensions
  DATA   = load('DATA.mat','TARIFF');
  rawT   = DATA.TARIFF;                 % [S·N × N]
[Nrows, N] = size(rawT);
S = Nrows / N;  assert(mod(S,1)==0,'mylambdaj: bad TARIFF shape.');
TARIFFs = permute(reshape(rawT',[N,N,S]), [2,1,3]);

  %% 2) Set bounds and prohibitive flags
  isProhib = (targetTariff_j > 2.25);        % identify prohibitive targets
  UB       = max(targetTariff_j + 0.03, 2.25);% upper bound for MFN tariff

  %% 3) Initialize guesses and storage
  guessMat   = [targetTariff_j, (1:S)', ones(S,1)];
  mfnGuess   = 0.1 * ones(S,1);
  posFlag    = ones(S,1);
  maxIter    = 250;
  progress   = zeros(S, 7, maxIter);
  converged  = false;
  iter       = 0;

  %% 4) Iteratively solve for lambda_j
  while ~converged && iter < maxIter
    % a) Extract current guess for lambda_j
    Lmat          = ones(N, S);
    Lmat(j, :)    = guessMat(:,3)';

    % b) Compute MFN-optimal tariffs for given lambda
    [mfnTariff_j, govWelfare] = mymfnoptimaltariffj(j, TARIFFs, Lmat, 0, UB, mfnGuess);

    % --- c) Measure deviation and build CHECK for logging ------------------
    tol = 1e-6;
    CHECK = sortrows([targetTariff_j, (1:S)', mfnTariff_j, isProhib, UB], [1 2]);
    
    % --- d/e) Update lambda guess toward targets (simple sign step) --------
    dev = mfnTariff_j - targetTariff_j;      % >0 means MFN > target
    dir = -sign(dev);                        % move to reduce gap
    
    % step size schedule: small, decaying
    base = 0.02; decay = 0.98;
    step = base * (decay^iter);
    
    KEEP = guessMat(:,3);                    % previous λ_j,s
    lambda_new = KEEP + step .* dir;         % step in direction
    lambda_new = max(0, lambda_new);         % enforce non-negativity
    lambda_new = lambda_new / max(eps, mean(lambda_new));  % normalise to mean 1
    
    CORRECTION = lambda_new - KEEP;
    guessMat(:,3) = lambda_new;              % apply update
    
    % --- f) Record progress snapshot (S×7 = S×5 + S×1 + S×1) --------------
    progress(:, :, iter+1) = [CHECK, CORRECTION, KEEP];

    % g) Update mfnGuess for next iteration
    if iter < 10
        mfnGuess = mfnTariff_j;                  % use the last solved MFN tariff
    else
        mfnGuess = targetTariff_j;               % fallback: target (S×1)
        if j == 5                                 % keep your Japan special case if wanted
            mfnGuess = targetTariff_j;
        end
    end

    % h) Increment iteration counter
    iter = iter + 1;
  end

  %% 5) Final selection if not fully converged
  if ~converged
    % pick the iterate with minimum RSS deviation
    rssVec = squeeze(sum(progress(:,5,1:iter).^2))';
    [~, bestIdx] = min(rssVec);
    bestSlice    = progress(:,:,bestIdx);
    mfnTariff_j  = bestSlice(:,3);
    lambda_j     = bestSlice(:,7);
  else
    lambda_j = guessMat(:,3);
  end
end
