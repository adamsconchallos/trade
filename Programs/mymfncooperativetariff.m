function [MFNCOOPERATIVETARIFFs, GOVERNMENTWELFAREHAT, EXPENDITUREHAT, WAGEHAT, ceq] = ...
    mymfncooperativetariff(LAMBDA, LBj, UBj, LIBONLY, MFNCOOPERATIVETARIFFGUESSs)
% Compute cooperative MFN tariffs (dimension-agnostic in N, S).
% Inputs:
%   LAMBDA  : N×S political weights (λ). Baseline: ones(N,S).
%   LBj,UBj : scalar bounds in levels (not logs) for MFN import tariffs.
%   LIBONLY : 0 allow ↑/↓, 1 only reductions (upper bound forced to 1 in transformed space).
%   MFNCOOPERATIVETARIFFGUESSs (optional): N×N×S bilateral guess (diag=0). 
% Outputs:
%   MFNCOOPERATIVETARIFFs: N×N×S, zeros on diagonal.
%   GOVERNMENTWELFAREHAT, EXPENDITUREHAT, WAGEHAT: N×1.
%   ceq: constraint residual vector at optimum.

  % Globals used by the shared objective/constraint
  global ARG SCALING TAU
  mycalculations;                      % defines N,S,TARIFFs,SIGMA, etc.

  % Basic checks
  assert(isscalar(LBj) && isscalar(UBj) && isscalar(LIBONLY),'LBj/UBj/LIBONLY must be scalars.');
  assert(LIBONLY==0 || LIBONLY==1,'LIBONLY must be 0 or 1.');
  [N1,N2,S] = size(TARIFFs);  assert(N1==N2,'TARIFFs must be N×N×S');  N = N1;

  ARG     = LAMBDA - 1;       % pass λ in levels via ARG+1 downstream
  SCALING = 0.1;

  % ---------- 1) Build MFN import-tariff guess (S×N) ----------
  if nargin>4 && ~isempty(MFNCOOPERATIVETARIFFGUESSs)
      % coerce to N×N×S if user passed a different but compatible shape
      if ~isequal(size(MFNCOOPERATIVETARIFFGUESSs), [N N S])
          assert(numel(MFNCOOPERATIVETARIFFGUESSs)==N*N*S, ...
                 'Guess has %d elements; expected %d for N×N×S.', ...
                  numel(MFNCOOPERATIVETARIFFGUESSs), N*N*S);
          MFNCOOPERATIVETARIFFGUESSs = reshape(MFNCOOPERATIVETARIFFGUESSs, [N N S]);
      end
      % MFN import guess = average across sources i≠j for each importer j and sector s
      MFNCOOPERATIVEIMPTARIFFGUESS = zeros(S, N);
      for j = 1:N
          MFNCOOPERATIVEIMPTARIFFGUESS(:, j) = reshape( ...
              mean(MFNCOOPERATIVETARIFFGUESSs([1:j-1, j+1:N], j, :), 1), [S,1]);
      end
  else
      if isequal(LAMBDA, ones(N,S))
          INVMARKUP = (SIGMA - 1) ./ SIGMA;              % S×1
          TEMP = 0.25 * (INVMARKUP ./ mean(INVMARKUP) - 1);
          MFNCOOPERATIVEIMPTARIFFGUESS = repmat(TEMP, 1, N);   % S×N
          % Build a bilateral cube guess with zero diagonals
          MFNCOOPERATIVETARIFFGUESSs = zeros(N,N,S);
          for j=1:N
              col = reshape(TEMP,[1 1 S]);                       % 1×1×S
              MFNCOOPERATIVETARIFFGUESSs(:,j,:) = repmat(col,[N,1,1]);
              MFNCOOPERATIVETARIFFGUESSs(j,j,:) = 0;
          end
      else
          MFNCOOPERATIVEIMPTARIFFGUESS = zeros(S, N);
          MFNCOOPERATIVETARIFFGUESSs   = zeros(N, N, S);
      end
  end
  MFNCOOPERATIVEIMPTARIFFGUESS = reshape(MFNCOOPERATIVEIMPTARIFFGUESS, S, N);  % S×N

  % ---------- 2) Transformations for fmincon ----------
  % MFN factual import tariffs (S×N): average over partners i≠j
  MFNIMPTARIFF = zeros(S,N);
  for j=1:N
      MFNIMPTARIFF(:,j) = reshape(mean(TARIFFs([1:j-1, j+1:N], j, :), 1), [S,1]);
  end

  TAU = 1 + SCALING * MFNIMPTARIFF;                      % S×N
  CGUESSTRANSF = (1 + SCALING * MFNCOOPERATIVEIMPTARIFFGUESS) ./ TAU;  % S×N
  LBCTRANSF    = (1 + SCALING * LBj) ./ TAU;                              
  if LIBONLY==0
      UBCTRANSF = (1 + SCALING * UBj) ./ TAU;
  else
      UBCTRANSF = ones(size(TAU));  % only reductions
  end
  CGUESSTRANSF = CGUESSTRANSF(:);   % (N*S)×1
  LBCTRANSF    = LBCTRANSF(:);
  UBCTRANSF    = UBCTRANSF(:);

  % ---------- 3) Initial GE point ----------
  [G0, ~, W0, ~, ~, X0] = mycounterfactuals(MFNCOOPERATIVETARIFFGUESSs, zeros(N,1), LAMBDA);
  G0 = G0(:);  X0 = X0(:);  W0 = W0(:);
  % stack decision vector x = [G; X; W; C_transf]
  GUESS = [G0; X0; W0; CGUESSTRANSF];
  % bounds
  LB = [0.1*ones(N,1); 0.1*ones(2*N,1); LBCTRANSF];
  UB = [10*ones(N,1);  10*ones(2*N,1);  UBCTRANSF];

  % ---------- 4) Solve ----------
  tol  = 1e-10;
  opts = optimset('algorithm','active-set','display','off', ...
                  'TolFun',tol,'TolX',tol,'TolCon',tol, ...
                  'MaxFunEvals',inf,'MaxIter',5000);
  x = fmincon(@myfuncoop, GUESS, [], [], [], [], LB, UB, @mymfnconcoop, opts);
  [~, ceq] = mymfnconcoop(x);   % residuals at optimum

  % ---------- 5) Unpack ----------
  GOVERNMENTWELFAREHAT = x(1:N);
  EXPENDITUREHAT       = x(N+1:2*N);
  WAGEHAT              = x(2*N+1:3*N);
  MFNCOOP_TRANSF       = x(3*N+1:end);                       % (N*S)×1

  % back-transform to MFN import tariffs S×N
  TEMP_import = (reshape(MFNCOOP_TRANSF, S, N) .* TAU - 1) ./ SCALING;  % S×N

  % inflate to bilateral N×N×S with zeros on diag
  MFNCOOPERATIVETARIFFs = zeros(N,N,S);
  for j=1:N
      col = reshape(TEMP_import(:,j), [1 1 S]);              % 1×1×S
      MFNCOOPERATIVETARIFFs(:,j,:) = repmat(col,[N,1,1]);
      MFNCOOPERATIVETARIFFs(j,j,:) = 0;
  end

  % ---------- 6) Save (legacy side-effects) ----------
  isBAS = isequal(LAMBDA, ones(N,S));
  if isBAS
      if LBj == 0
          RESTRICTEDMFNCOOPERATIVETARIFFBASs = MFNCOOPERATIVETARIFFs; %#ok<NASGU>
          save('RESTRICTEDMFNCOOPERATIVETARIFFBASs.mat','RESTRICTEDMFNCOOPERATIVETARIFFBASs');
      elseif LBj < 0
          UNRESTRICTEDMFNCOOPERATIVETARIFFBASs = MFNCOOPERATIVETARIFFs; %#ok<NASGU>
          save('UNRESTRICTEDMFNCOOPERATIVETARIFFBASs.mat','UNRESTRICTEDMFNCOOPERATIVETARIFFBASs');
      end
  else
      if LBj == 0
          RESTRICTEDMFNCOOPERATIVETARIFFPOLs = MFNCOOPERATIVETARIFFs; %#ok<NASGU>
          save('RESTRICTEDMFNCOOPERATIVETARIFFPOLs.mat','RESTRICTEDMFNCOOPERATIVETARIFFPOLs');
      elseif LBj < 0
          UNRESTRICTEDMFNCOOPERATIVETARIFFPOLs = MFNCOOPERATIVETARIFFs; %#ok<NASGU>
          save('UNRESTRICTEDMFNCOOPERATIVETARIFFPOLs.mat','UNRESTRICTEDMFNCOOPERATIVETARIFFPOLs');
      end
  end
end
