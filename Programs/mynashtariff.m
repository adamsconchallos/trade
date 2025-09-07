function NASHTARIFFs = mynashtariff(LAMBDA, LBj, UBj, NASHTARIFFGUESSs)
% MYNASHTARIFF  Compute bilateral Nash tariffs (non-MFN case).
%
%   NASHTARIFFs = mynashtariff(LAMBDA, LBj, UBj, NASHTARIFFGUESSs)
%
% Inputs
%   LAMBDA               [N×S]    political-economy weights
%   LBj, UBj             scalar (or [S×1]) lower/upper bounds per importer
%   NASHTARIFFGUESSs     [N×N×S]  initial bilateral tariff cube
%                                 (column j is importer j’s tariffs to all i)
%
% Output
%   NASHTARIFFs          [N×N×S]  fixed point of best responses
%
% Notes
%   For each importer j, myoptimaltariffj returns a vector of length (N-1)*S
%   containing the off-diagonal bilateral tariffs by sector. We embed that
%   vector back as column j of the [N×N×S] cube with zeros on the diagonal.
%
% Parallelism
%   Uses an existing parallel pool if one is open (parpool / gcp).
%   Otherwise runs serially. No matlabpool calls.

  % ---- Validate shapes -------------------------------------------------
  [N1, N2, S] = size(NASHTARIFFGUESSs);
  assert(N1 == N2, 'mynashtariff: guess cube must be N×N×S.');
  N = N1;
  assert(isequal(size(LAMBDA), [N, S]), 'mynashtariff: LAMBDA must be N×S.');

  % ---- Initialise state ------------------------------------------------
  TARIFF_curr = NASHTARIFFGUESSs;              % current cube  [N×N×S]

  % Build the (N-1)S×N initial off-diagonal guesses for each importer j
  OPTIMALTARIFFGUESS = zeros((N-1)*S, N);
  for j = 1:N
      tmp = TARIFF_curr(:, j, :);              % [N×1×S]
      tmp(j,:,:) = [];                         % drop diagonal row i=j
      OPTIMALTARIFFGUESS(:, j) = reshape(tmp, (N-1)*S, 1); % flatten to vector for optimiser
  end

  % ---- Iteration controls ---------------------------------------------
  maxit = 24;
  tol   = 1e-2;
  err   = inf;
  t     = 0;

  % ---- Iterate best responses -----------------------------------------
  while (t < maxit) && (err > tol)
      TARIFF_next = TARIFF_curr;

      % Use parallel if a pool exists; otherwise serial
      if ~isempty(gcp('nocreate'))
          parfor j = 1:N
              colj = best_response_column(j, TARIFF_curr, LAMBDA, LBj, UBj, OPTIMALTARIFFGUESS(:, j), S, N);
              TARIFF_next(:, j, :) = colj;
          end
      else
          for j = 1:N
              colj = best_response_column(j, TARIFF_curr, LAMBDA, LBj, UBj, OPTIMALTARIFFGUESS(:, j), S, N);
              TARIFF_next(:, j, :) = colj;
          end
      end

      % Convergence: mean absolute change over the cube
      err = mean(abs(TARIFF_next(:) - TARIFF_curr(:)));

      % Advance
      TARIFF_curr = TARIFF_next;
      t = t + 1;
  end

  NASHTARIFFs = TARIFF_curr;
end


% -------------------- local subfunction --------------------------------
function colj = best_response_column(j, cube, LAMBDA, LBj, UBj, guess_j, S, N)
% Compute importer j’s best-response bilateral tariffs and embed as column j.

  % myoptimaltariffj returns ((N-1)*S)×1 vector (off-diagonal bilaterals)
  opt_j = myoptimaltariffj(j, cube, LAMBDA, LBj, UBj, guess_j);

  % Embed into column j of the [N×N×S] cube with zero diagonal
  tmp  = reshape(opt_j, [N-1, 1, S]);          % (N-1)×1×S
  colj = [tmp(1:j-1, :, :) ; zeros(1,1,S) ; tmp(j:end, :, :)];  % N×1×S
end
