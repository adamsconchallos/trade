function MFNNASHTARIFFs = mymfnnashtariff(LAMBDA, LBj, UBj, MFNNASHTARIFFGUESSs)
% MYMFNNASHTARIFF  Nash equilibrium of MFN tariffs (per importer j)
%
%   MFNNASHTARIFFs = mymfnnashtariff(LAMBDA, LBj, UBj, MFNNASHTARIFFGUESSs)
%
% Inputs
%   LAMBDA                 [N×S]  political-economy weights
%   LBj, UBj               scalars (or [S×1]) lower/upper bounds for tariffs
%   MFNNASHTARIFFGUESSs    [N×N×S] initial cube, column j is importer j’s MFN column
%
% Output
%   MFNNASHTARIFFs         [N×N×S] fixed point of MFN best responses
%
% Notes
%   For importer j: the MFN schedule is a [S×1] vector applied to all partners i≠j;
%   the bilateral column j in the cube is that vector replicated with zeros on the diagonal.

  % ---- Validate shapes -------------------------------------------------
  [N1, N2, S] = size(MFNNASHTARIFFGUESSs);
  assert(N1==N2, 'mymfnnashtariff: guess cube must be N×N×S.');
  N = N1;
  assert(isequal(size(LAMBDA), [N, S]), 'mymfnnashtariff: LAMBDA must be N×S.');

  % ---- Initialise state ------------------------------------------------
  TARIFFOTHER_curr = MFNNASHTARIFFGUESSs;   % current cube [N×N×S]

  % MFN import tariff guess for each importer j: mean over i≠j of column j
  MFNIMP_curr = zeros(S, N);
  for j = 1:N
      tmp = TARIFFOTHER_curr(:, j, :);      % [N×1×S]
      tmp(j,:,:) = [];                      % drop diagonal
      MFNIMP_curr(:, j) = reshape(mean(tmp, 1), [S, 1]);
  end

  % ---- Iteration controls ---------------------------------------------
  maxit = 24;
  tol   = 1e-2;
  err   = inf;
  t     = 0;

  % ---- Iterate best responses -----------------------------------------
  while (t < maxit) && (err > tol)
      TARIFFOTHER_next = TARIFFOTHER_curr;
      MFNIMP_next      = MFNIMP_curr;

      % Use parallel if a pool exists; otherwise fall back to serial loop
      if ~isempty(gcp('nocreate'))
          parfor j = 1:N
              [colj, mfnj] = update_column(j, TARIFFOTHER_curr, LAMBDA, LBj, UBj, MFNIMP_curr(:, j), S, N);
              TARIFFOTHER_next(:, j, :) = colj;
              MFNIMP_next(:, j)         = mfnj;
          end
      else
          for j = 1:N
              [colj, mfnj] = update_column(j, TARIFFOTHER_curr, LAMBDA, LBj, UBj, MFNIMP_curr(:, j), S, N);
              TARIFFOTHER_next(:, j, :) = colj;
              MFNIMP_next(:, j)         = mfnj;
          end
      end

      % Convergence check (mean absolute change over the cube)
      err = mean(abs(TARIFFOTHER_next(:) - TARIFFOTHER_curr(:)));

      % Advance
      TARIFFOTHER_curr = TARIFFOTHER_next;
      MFNIMP_curr      = MFNIMP_next;
      t = t + 1;
  end

  MFNNASHTARIFFs = TARIFFOTHER_curr;
end


% -------------------- local subfunction --------------------------------
function [colj, mfn_j] = update_column(j, cube, LAMBDA, LBj, UBj, mfnGuess_j, S, N)
% Compute importer j’s MFN best response and embed as column j of the cube.

  % Best-response MFN tariff for importer j (S×1)
  mfn_j = mymfnoptimaltariffj(j, cube, LAMBDA, LBj, UBj, mfnGuess_j);

  % Replicate MFN vector down the column and zero the diagonal
  tmp  = repmat(reshape(mfn_j, [1, 1, S]), [N-1, 1, 1]);  % (N-1)×1×S
  colj = [tmp(1:j-1, :, :) ; zeros(1,1,S) ; tmp(j:end, :, :)];  % N×1×S
end
