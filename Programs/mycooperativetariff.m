%------------------------------------------
% Compute cooperative (non-MFN) tariffs
%------------------------------------------
function [COOPERATIVETARIFFs, GOVERNMENTWELFAREHAT, EXPENDITUREHAT, WAGEHAT, ceq] = ...
    mycooperativetariff(LAMBDA, LBj, UBj, LIBONLY, COOPERATIVETARIFFGUESSs)

    % Globals read by the objective/constraint (legacy design)
    global ARG SCALING TAU

    % --- Preliminaries & data ---
    mycalculations;                        % must define N,S,TARIFFs,SIGMA,…
    ARG = LAMBDA - 1;                      % political-economy weights in levels → deviations

    % --- Input checks (fail fast) ---
    assert(isscalar(LBj) && isscalar(UBj) && isscalar(LIBONLY), 'LBj, UBj, LIBONLY must be scalars.');
    assert(LIBONLY==0 || LIBONLY==1, 'LIBONLY must be 0 or 1.');
    assert(isequal(size(LAMBDA), [N, S]), 'LAMBDA must be N×S.');
    isBAS = isequal(LAMBDA, ones(N,S));    % BAS=all ones; POL otherwise

    % -----------------------------------------------------------
    % 1) Build initial guess cube (N×N×S) and off-diagonal stack
    % -----------------------------------------------------------
    if nargin > 4 && ~isempty(COOPERATIVETARIFFGUESSs)
        % Use user-supplied bilateral guess; coerce to N×N×S
        COOPERATIVETARIFFGUESSs = alignToNNScube(COOPERATIVETARIFFGUESSs, TARIFFs);
    else
        % No guess supplied → replicate Ossa’s baseline logic
        if isBAS
            % Target elimination of markup distortions (σ/(σ−1) markups)
            INVMARKUP = (SIGMA - 1) ./ SIGMA;              % S×1
            TEMP = 0.25 * (INVMARKUP ./ mean(INVMARKUP) - 1);  % S×1
            % Build bilateral guess with zero diagonals
            COOPERATIVETARIFFGUESSs = zeros(N, N, S);
            base = repmat(reshape(TEMP, [1 1 S]), [N-1, 1, 1]); % (N−1)×1×S
            for j = 1:N
                COOPERATIVETARIFFGUESSs(:, j, :) = [base(1:j-1,1,:); zeros(1,1,S); base(j:end,1,:)];
            end
        else
            COOPERATIVETARIFFGUESSs = zeros(N, N, S);
        end
    end

    % Extract importer-j off-diagonal import tariffs into a compact stack
    % COOPERATIVEIMPTARIFFGUESSs has size (N−1)×N×S (rows = origins i≠j)
    COOPERATIVEIMPTARIFFGUESSs = zeros(N-1, N, S);
    for j = 1:N
        COOPERATIVEIMPTARIFFGUESSs(:, j, :) = COOPERATIVETARIFFGUESSs([1:j-1, j+1:N], j, :);
    end

    % ----------------------------------------------------------------
    % 2) Build initial GE state for fmincon decision vector [G;X;W;C]
    % ----------------------------------------------------------------
    [G0, ~, W0, ~, ~, X0] = mycounterfactuals(COOPERATIVETARIFFGUESSs, zeros(N,1), LAMBDA);

    % Enforce column shapes to avoid vertical-concat errors
    G0 = G0(:);  X0 = X0(:);  W0 = W0(:);
    assert(numel(G0)==N && numel(X0)==N && numel(W0)==N, ...
        'mycounterfactuals returned unexpected sizes: G=%s, X=%s, W=%s.', ...
        mat2str(size(G0)), mat2str(size(X0)), mat2str(size(W0)));

    % ----------------------------------------------------------------
    % 3) Transformation for tariff variables (conditioning trick)
    % ----------------------------------------------------------------
    SCALING = 0.1;

    % Stack factual import tariffs by importer j and origins i≠j
    IMPTARIFFs = zeros(N-1, N, S);
    for j = 1:N
        IMPTARIFFs(:, j, :) = TARIFFs([1:j-1, j+1:N], j, :);
    end
    TAU = 1 + SCALING * IMPTARIFFs;   % same size as IMPTARIFFs

    % Transform guesses and bounds to work on tariff *changes* scaled by TAU
    CGUESSTRANSF = ((1 + SCALING * COOPERATIVEIMPTARIFFGUESSs) ./ TAU);
    LBCTRANSF    = ((1 + SCALING * LBj) ./ TAU);
    if LIBONLY == 0
        UBCTRANSF = ((1 + SCALING * UBj) ./ TAU);
    else
        UBCTRANSF = ones(size(TAU));  % only reductions allowed
    end

    % Vectorize to columns: [(N−1)·N·S]×1
CGUESSTRANSF = CGUESSTRANSF(:);
LBCTRANSF    = LBCTRANSF(:);
UBCTRANSF    = UBCTRANSF(:);

    % Decision vector and simple bounds:
    % x = [GOVERNMENTWELFAREHAT (N); EXPENDITUREHAT (N); WAGEHAT (N); C_transf ((N−1)·N·S)]
    GUESS = [G0; X0; W0; CGUESSTRANSF];
    LB    = [0.1*ones(N,1); 0.1*ones(2*N,1); LBCTRANSF];
    UB    = [10*ones(N,1);  10*ones(2*N,1);  UBCTRANSF];

    % -----------------------------
    % 4) Solve cooperative program
    % -----------------------------
    tol = 1e-10;
    opts = optimoptions('fmincon', ...
        'Algorithm','interior-point', ...           % or 'sqp' if preferred
        'Display','off', ...
        'OptimalityTolerance',tol, ...
        'StepTolerance',tol, ...
        'ConstraintTolerance',tol, ...
        'MaxFunctionEvaluations',1e6, ...
        'MaxIterations',5000);

    % Legacy line if you insist on the old API:
    % opts = optimset('algorithm','active-set','display','off','TolFun',tol,'TolX',tol,'TolCon',tol,'MaxFunEvals',inf,'MaxIter',5000);
    x = fmincon(@myfuncoop, GUESS, [], [], [], [], LB, UB, @myconcoop, opts);

    % Constraint residuals (diagnostic)
    [~, ceq] = myconcoop(x);

    % -----------------------------
    % 5) Unpack decision variables
    % -----------------------------
    GOVERNMENTWELFAREHAT = x(1:N);
    EXPENDITUREHAT       = x(N+1:2*N);
    WAGEHAT              = x(2*N+1:3*N);
    COOP_TRANSF          = x(3*N+1:end);                % transformed tariffs

    % Back-transform to levels: off-diagonal stack → bilateral cube with zeros on diagonal
    TEMP = (reshape(COOP_TRANSF, [N-1, N, S]) .* TAU - 1) ./ SCALING;   % (N−1)×N×S
    COOPERATIVETARIFFs = zeros(N, N, S);
    for j = 1:N
        COOPERATIVETARIFFs(:, j, :) = [TEMP(1:j-1, j, :); zeros(1,1,S); TEMP(j:end, j, :)];
    end

    % -----------------------------
    % 6) Save side-effect files
    % -----------------------------
    if isBAS
        if LBj == 0
            RESTRICTEDCOOPERATIVETARIFFBASs = COOPERATIVETARIFFs; %#ok<NASGU>
            save('RESTRICTEDCOOPERATIVETARIFFBASs.mat','RESTRICTEDCOOPERATIVETARIFFBASs');
        elseif LBj < 0
            UNRESTRICTEDCOOPERATIVETARIFFBASs = COOPERATIVETARIFFs; %#ok<NASGU>
            save('UNRESTRICTEDCOOPERATIVETARIIFFBASs.mat','UNRESTRICTEDCOOPERATIVETARIFFBASs');  % file name can match your pipeline
        end
    else
        if LBj == 0
            RESTRICTEDCOOPERATIVETARIFFPOLs = COOPERATIVETARIFFs; %#ok<NASGU>
            save('RESTRICTEDCOOPERATIVETARIFFPOLs.mat','RESTRICTEDCOOPERATIVETARIFFPOLs');
        elseif LBj < 0
            UNRESTRICTEDCOOPERATIVETARIFFPOLs = COOPERATIVETARIFFs; %#ok<NASGU>
            save('UNRESTRICTEDCOOPERATIVETARIFFPOLs.mat','UNRESTRICTEDCOOPERATIVETARIFFPOLs');
        end
    end
end

% ---------- Local utility ----------
function G = alignToNNScube(X, TARIFFs)
    tgt = size(TARIFFs);                      % [N N S]
    if isequal(size(X), tgt), G = X; return; end
    if numel(X) ~= prod(tgt)
        error('Guess has %d elements; expected %d for N×N×S.', numel(X), prod(tgt));
    end
    G = reshape(X, tgt);
end
