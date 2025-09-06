%=======================================================================
% results.m  –  Master driver to replicate numerical results in
%               Ossa (2014) “Trade Wars and Trade Talks with Data”
%               *American Economic Review* 104(12): 4104–46.
%-------------------------------------------------------------------------
% SCRIPT ↔ PAPER MAPPING
% ┌──────────┬───────────────────────────────────────────────────────────┐
% │  Section │ MATLAB Block                                              │
% ├──────────┼───────────────────────────────────────────────────────────┤
% │ §III     │ PART 1: Data work & calibration (eqs. 6–9 in levels)      │
% │ §IV.B    │ PART 1E: Political‑economy λ calibration                  │
% │ §V       │ PART 2: Unilateral optimal tariffs (Table 3)              │
% │ §VI      │ PART 3: Trade Wars – Nash equilibrium (Table 4)           │
% │ §VII     │ PART 4: Trade Talks – efficient cooperation (Fig.2, Tbl5) │
% │ §VIII    │ PART 5: Robustness / misc. (Online App.)                  │
% └──────────┴───────────────────────────────────────────────────────────┘
% Note: All GE hat-systems (eqs. 10–13) run via:
%    mycounterfactuals() → mysolvedconditions() → myconditions().
%=======================================================================

%% --- USER SETTINGS ----------------------------------------------------
projectRoot = "C:\Users\adams\OneDrive - University of Florida\Research\Emran_trade\Data_MS_AER_2012_0527";
versions    = {"Main"};
useParallel = true;
useToyData = true;            % ← flip to false for the full dataset
dataSub    = ternary(useToyData, 'toy CSV', 'CSV files');
chs = @(varargin) cd(fullfile(projectRoot,varargin{:}));

%% --- SET PATH & POOL --------------------------------------------------
addpath(projectRoot);
addpath(fullfile(projectRoot,'Programs'));
if useParallel
    delete(gcp('nocreate'));
    parpool;
end

%% ---------------------------------------------------------------------
% PART 1 – §III: Data & calibration (eqs. 6–9)
%-----------------------------------------------------------------------
chs(); clearvars -except projectRoot versions useParallel chs dataSub useToyData;
for v = 1:numel(versions)
    vers = versions{v}; 
    dataDir = fullfile(projectRoot,'Data',vers);
    mkdir_if(dataDir);          % ← create Data/<vers> if needed
    chs('Data',vers);
    copyfile(fullfile(projectRoot,'Data',dataSub,'*.csv'), dataDir, 'f');
    read   = @(f) csvread(fullfile(projectRoot,'Data',dataSub,f));
    TRADE        = read('trade_ijs.csv');
    TARIFF       = read('tariff_ijs.csv');
    SIGMAALT     = read('sigma_s.csv');
    SIGMA        = SIGMAALT(:,v);
    TARGETTARIFF = read('tariff_noncoop.csv');
    LABELS       = strtrim(importdata(fullfile(projectRoot,'Data',dataSub,'labels.csv')));
    save DATA.mat TRADE TARIFF SIGMA TARGETTARIFF LABELS SIGMAALT;

    % ← START USING cfg.dims HERE
    [N,S] = cfg.dims.loadDims(projectRoot, dataSub);
    % Baseline calibration (eqs. 6–9)
    mycalculations;
    % Factual GE‐consistent trade (eqs. 10–13)
    [~,~,~,TRADECs] = mycounterfactuals(TARIFFs, zeros(N,1), ones(N,S));
    TRADE = reshape(permute(TRADECs, [1 3 2]), [N*S, N, 1]);
    save DATA.mat -append TRADE;

    % Political‐economy λ (§ IV.B)
    mycalculations;
    for jj = [1,4,6]
        TARGETTARIFF(:, jj) = MEANTARIFF(:, jj);
    end
    LAMBDAPOL   = zeros(N, S);
    MFNOPTIMALTARIFFPOL = zeros(S, N);
    parfor j = 1:N
        [Lj, Mf] = mylambdaj(j, TARGETTARIFF(:, j));
        LAMBDAPOL(j, :) = Lj';
        MFNOPTIMALTARIFFPOL(:, j) = Mf;
    end
    LAMBDABAS = ones(N, S);
    save DATA.mat -append LAMBDAPOL LAMBDABAS MFNOPTIMALTARIFFPOL;
end

info = whos('-file','DATA.mat');
names = {info.name};
assert( all(ismember({'LAMBDABAS','LAMBDAPOL'}, names)), ...
        'DATA.mat is missing LAMBDABAS/LAMBDAPOL – PART 1 did not complete.');
%% ---------------------------------------------------------------------
% PART 2 – §V: Unilateral optimal tariffs (Table 3)  ~85m
% ----------------------------------------------------------------------
for v = 1:numel(versions)
    vers     = versions{v};
    dataPath = fullfile(projectRoot, 'Data', vers, 'DATA.mat');
    tmp      = load(dataPath, 'TARIFF', 'LAMBDABAS', 'LAMBDAPOL');
    TARIFF   = tmp.TARIFF;
    LAMBDABAS = tmp.LAMBDABAS;
    LAMBDAPOL = tmp.LAMBDAPOL;

    N = size(TARIFF, 2);
    S = size(TARIFF, 1) / N;
    TARIFFs = permute(reshape(TARIFF', [N, N, S]), [2, 1, 3]);

    MFNOPTIMALTARIFFBAS = zeros(S, N);
    MFNOPTIMALTARIFFPOL = zeros(S, N);
    OPTIMALTARIFFBAS    = zeros((N-1)*S, N);
    OPTIMALTARIFFPOL    = zeros((N-1)*S, N);

    parfor j = 1:N
        MFNOPTIMALTARIFFBAS(:, j) = mymfnoptimaltariffj(j, TARIFFs, LAMBDABAS, 0, 10);
        MFNOPTIMALTARIFFPOL(:, j) = mymfnoptimaltariffj(j, TARIFFs, LAMBDAPOL, 0, 10);

        OPTIMALTARIFFBAS(:, j) = myoptimaltariffj( ...
            j, TARIFFs, LAMBDABAS, 0, 10, ...
            reshape(repmat(MFNOPTIMALTARIFFBAS(:, j)', N-1, 1), (N-1)*S, 1) ...
        );
        OPTIMALTARIFFPOL(:, j) = myoptimaltariffj( ...
            j, TARIFFs, LAMBDAPOL, 0, 10, ...
            reshape(repmat(MFNOPTIMALTARIFFPOL(:, j)', N-1, 1), (N-1)*S, 1) ...
        );
    end

    resF = fullfile(projectRoot, 'Results', 'Optimal tariffs', vers);
    if ~exist(resF, 'dir'), mkdir(resF); end
    save(fullfile(resF, 'MFNOPTIMALTARIFFBAS.mat'), 'MFNOPTIMALTARIFFBAS');
    save(fullfile(resF, 'MFNOPTIMALTARIFFPOL.mat'), 'MFNOPTIMALTARIFFPOL');
    save(fullfile(resF, 'OPTIMALTARIFFBAS.mat'), 'OPTIMALTARIFFBAS');
    save(fullfile(resF, 'OPTIMALTARIFFPOL.mat'), 'OPTIMALTARIFFPOL');
end
%% ---------------------------------------------------------------------
% PART 3 – §VI: Trade Wars (Nash eq., Table 4)  ~96m
% ----------------------------------------------------------------------
clearvars -except projectRoot versions useParallel chs dataSub useToyData; 

for v = 1:numel(versions)
    vers = versions{v};  
    resV=fullfile(projectRoot,'Results','Trade wars',vers);
    if ~exist(resV,'dir'), mkdir(resV); end
    copyfile(fullfile(projectRoot,'Data',vers,'DATA.mat'),resV);
    copyfile(fullfile(projectRoot,'Results','Optimal tariffs',vers,'*.mat'),resV);
    copyfile(fullfile(projectRoot,'Programs','*.m'),resV);
    chs('Results','Trade wars',vers);
    mycalculations;

    load MFNOPTIMALTARIFFBAS.mat
    MFNNASHTARIFFBASs= mymfnnashtariff(LAMBDABAS,0,10,restack(MFNOPTIMALTARIFFBAS)); save MFNNASHTARIFFBASs.mat MFNNASHTARIFFBASs
    load MFNOPTIMALTARIFFPOL.mat
    MFNNASHTARIFFPOLs= mymfnnashtariff(LAMBDAPOL,0,10,restack(MFNOPTIMALTARIFFPOL));save MFNNASHTARIFFPOLs.mat MFNNASHTARIFFPOLs

    NASHTARIFFBASs= mynashtariff(LAMBDABAS,0,10,MFNNASHTARIFFBASs);
    NASHTARIFFPOLs= mynashtariff(LAMBDAPOL,0,10,MFNNASHTARIFFPOLs);
    save NASHTARIFFBASs.mat NASHTARIFFBASs
    save NASHTARIFFPOLs.mat NASHTARIFFPOLs

    %delete('*.m','*TARIFF*.mat','DATA.mat'); chs();
end

%% ---------------------------------------------------------------------
% PART 4 – §VII: Trade Talks (Efficient coop., Fig2/Tbl5)
% ----------------------------------------------------------------------
clearvars -except projectRoot versions useParallel chs dataSub useToyData;

for v = 1:numel(versions)
    vers = versions{v};

    % 0) Prepare scenario folders and copy code/base data (no writes by workers)
    prepTradeTalkFolders(projectRoot, vers);

    % 1) SEED SCENARIOS (sequential; each writes DATA.mat exactly once)
    seedScenario(projectRoot, vers, 'Nash_Bas');  % DATA.mat = Nash_Bas baseline
    seedScenario(projectRoot, vers, 'Nash_Pol');  % DATA.mat = Nash_Pol baseline
    seedScenario(projectRoot, vers, 'Free');      % DATA.mat = zero tariffs
    seedScenario(projectRoot, vers, 'Fact');      % DATA.mat = factual

    % 2) COMPUTE (parallel; read-only DATA.mat)
    poolOpenedHere = false;
    if useParallel && isempty(gcp('nocreate')), parpool; poolOpenedHere = true; end

    parfor sc = 1:6
        runCoop_compute_only(projectRoot, vers, sc);
    end

    if poolOpenedHere, delete(gcp('nocreate')); end

    % 3) Optional post-process and clean
    if strcmp(vers,'Main')
        polishCoop(projectRoot, vers);
    end
    cleanMFiles(fullfile(projectRoot,'Results','Trade talks',vers));
end


%% ---------------------------------------------------------------------
% PART 5 – §VIII: Robustness / Online App.  ~contextual tweaks
% ----------------------------------------------------------------------
recalcRestrictedCoop(projectRoot,versionsFull);

%% ========= LOCAL FUNCTIONS ===========================================
function mkdir_if(p)
%MKDIR_IF  Create folder if needed. Add to MATLAB and worker paths when safe.
% Usage: mkdir_if(pathString)

    if ~exist(p,'dir'), mkdir(p); end

    % Avoid path shadowing: only add if the folder has no .m files
    hasMfiles = ~isempty(dir(fullfile(p,'*.m')));
    if ~hasMfiles
        % add on client if not already present
        if ~contains([path pathsep], [p pathsep], 'IgnoreCase', ispc)
            addpath(p);
        end
        % add on workers if a pool exists
        pool = gcp('nocreate');
        if ~isempty(pool)
            f = parfevalOnAll(@addpath,0,p);  % non-blocking broadcast
            wait(f);                           % ensure workers received it
        end
    end
end

function out = ternary(cond,a,b), if cond, out=a; else, out=b; end, end
function out = restack(vec)
% RESTACK
%   If vec is S×N (MFN by sector and importer), return N×N×S with
%   column j equal to MFN(:,j) for all i≠j and zeros on the diagonal.
%   If vec is (N-1)S×N (off-diagonal bilateral vector), return the
%   N×N×S cube by inserting zeros on the diagonal.

    [R, N] = size(vec);

    % Case A: MFN input S×N  → replicate across partners (zero diagonal)
    % Detect “MFN” by trying to interpret R as S (and not divisible by N-1)
    if mod(R, N-1) ~= 0
        S = R;                          % R is the number of sectors
        out = zeros(N, N, S);
        for j = 1:N
            mfn = reshape(vec(:,j), [1 1 S]);   % 1×1×S
            out(:, j, :) = repmat(mfn, [N 1 1]);% fill column j
            out(j, j, :) = 0;                   % zero diagonal
        end
        return
    end

    % Case B: bilateral vector input (N-1)S×N → insert zeros on diagonal
    S = R / (N - 1);
    assert(abs(S - round(S)) < eps, 'restack: incompatible first dimension.');

    out = zeros(N, N, S);
    for j = 1:N
        % take the (N-1)×1×S stack and insert a zero row at i=j
        tmp = reshape(vec(:, j), [N-1, 1, S]);
        out(:, j, :) = [tmp(1:j-1, :, :) ; zeros(1,1,S) ; tmp(j:end, :, :)];
    end
end

function prepTradeTalkFolders(root, vers)
% Create per-scenario folders and copy base DATA and code + needed Nash files.
    baseData = fullfile(root,'Data',vers,'DATA.mat');
    assert(exist(baseData,'file')==2, 'Missing base DATA: %s', baseData);

    subs = {'Nash_Bas','Nash_Pol','Free','Fact'};
    for k = 1:numel(subs)
        f = fullfile(root,'Results','Trade talks',vers,subs{k});
        if ~exist(f,'dir'), mkdir(f); end

        % Always bring pristine base DATA and the code
        copyfile(baseData, f);
        copyfile(fullfile(root,'Programs','*.m'), f);

        % For Nash_* also bring BOTH Nash tariff files from Trade wars
        if strcmp(subs{k}, 'Nash_Bas')
            srcs = { 'NASHTARIFFBASs.mat', 'MFNNASHTARIFFBASs.mat' };
        elseif strcmp(subs{k}, 'Nash_Pol')
            srcs = { 'NASHTARIFFPOLs.mat', 'MFNNASHTARIFFPOLs.mat' };
        else
            srcs = {};
        end
        for r = 1:numel(srcs)
            src = fullfile(root,'Results','Trade wars',vers,srcs{r});
            assert(exist(src,'file')==2, 'Missing file: %s', src);
            copyfile(src, f);
        end
    end
end


function seedScenario(root, vers, sub)
% Open target folder, build the intended T0, and write DATA.mat ONCE.
    folder = fullfile(root,'Results','Trade talks',vers,sub);
    cd(folder);

    % Load pristine base DATA from Data/<vers>/ to avoid stale local DATA
    base = load(fullfile(root,'Data',vers,'DATA.mat'));
    SIGMA     = base.SIGMA;
    TARIFF    = base.TARIFF;
    TRADE     = base.TRADE;
    LAMBDABAS = base.LAMBDABAS;
    LAMBDAPOL = base.LAMBDAPOL;

    % Save a fresh local DATA and derive N,S,TARIFFs
    save('DATA.mat','SIGMA','TARIFF','TRADE','LAMBDABAS','LAMBDAPOL','-mat');
    mycalculations;                 % defines N,S,TARIFFs from local DATA.mat

    switch sub
        case 'Nash_Bas'
            load NASHTARIFFBASs.mat
            T0 = NASHTARIFFBASs;
            [~,~,~,TRADECs] = mycounterfactuals(T0, zeros(N,1), LAMBDABAS);
            TARIFF = reshape(permute(T0,[1 3 2]), [N*S, N, 1]);
            TRADE  = reshape(permute(TRADECs,[1 3 2]), [N*S, N, 1]);
            save('DATA.mat','SIGMA','TARIFF','TRADE','LAMBDABAS','LAMBDAPOL','-mat');

        case 'Nash_Pol'
            load NASHTARIFFPOLs.mat
            T0 = NASHTARIFFPOLs;
            [~,~,~,TRADECs] = mycounterfactuals(T0, zeros(N,1), LAMBDAPOL);
            TARIFF = reshape(permute(T0,[1 3 2]), [N*S, N, 1]);
            TRADE  = reshape(permute(TRADECs,[1 3 2]), [N*S, N, 1]);
            save('DATA.mat','SIGMA','TARIFF','TRADE','LAMBDABAS','LAMBDAPOL','-mat');

        case 'Free'
            T0 = zeros(size(TARIFFs));      % N×N×S
            [~,~,~,TRADECs] = mycounterfactuals(T0, zeros(N,1), LAMBDABAS);
            TARIFF = reshape(permute(T0,[1 3 2]), [N*S, N, 1]);
            TRADE  = reshape(permute(TRADECs,[1 3 2]), [N*S, N, 1]);
            save('DATA.mat','SIGMA','TARIFF','TRADE','LAMBDABAS','LAMBDAPOL','-mat');

        case 'Fact'
            % Factual already saved; nothing else to write.
    end

    cd(root);
end

function runCoop_compute_only(root, vers, sc)
% Scenario mapping (as in the original script):
%   1: Nash_Bas  (BAS)  → unrestricted + restricted (BAS)
%   2: Fact      (BAS)  → unrestricted + restricted (BAS)
%   3: Free      (BAS)  → unrestricted + restricted (BAS)
%   4: Nash_Pol  (POL)  → restricted only (POL), guesses from Nash files
%   5: Fact      (POL)  → restricted only (POL), guess = TARIFFs
%   6: Free      (POL)  → restricted only (POL), guess = zeros

    subs   = {'Nash_Bas','Fact','Free','Nash_Pol','Fact','Free'};
    folder = fullfile(root,'Results','Trade talks',vers,subs{sc});
    cd(folder);

    mycalculations;  % builds N,S,TARIFFs and exposes LAMBDABAS/LAMBDAPOL

    if sc <= 3
        % BAS: unrestricted + restricted
        L  = LAMBDABAS;
        T0 = TARIFFs;

        UM = mymfncooperativetariff(L,-5,10,0,T0);
        U  = mycooperativetariff(L,-5,10,0,UM);

        RM = mymfncooperativetariff(L,0,10,0,max(UM,0));
        R  = mycooperativetariff(L,0,10,0,max(U,0));

        save('UNRESTRICTEDMFNCOOPERATIVETARIFFBASs.mat','UM');
        save('UNRESTRICTEDCOOPERATIVETARIFFBASs.mat','U');
        save('RESTRICTEDMFNCOOPERATIVETARIFFBASs.mat','RM');
        save('RESTRICTEDCOOPERATIVETARIFFBASs.mat','R');

    else
        % POL: restricted only
        L = LAMBDAPOL;

        if sc == 4
            % Nash_Pol: MFN guess = MFNNASHTARIFFPOLs, Coop guess = NASHTARIFFPOLs
            A  = load('MFNNASHTARIFFPOLs.mat','MFNNASHTARIFFPOLs');
            B  = load('NASHTARIFFPOLs.mat','NASHTARIFFPOLs');
            RM = mymfncooperativetariff(L,0,10,0,A.MFNNASHTARIFFPOLs);
            R  = mycooperativetariff(L,0,10,0,B.NASHTARIFFPOLs);

        elseif sc == 5
            % Fact (POL): guess = factual tariffs
            T0 = TARIFFs;
            RM = mymfncooperativetariff(L,0,10,0,T0);
            R  = mycooperativetariff(L,0,10,0,T0);

        else % sc == 6
            % Free (POL): guess = zeros
            T0 = zeros(size(TARIFFs));
            RM = mymfncooperativetariff(L,0,10,0,T0);
            R  = mycooperativetariff(L,0,10,0,T0);
        end

        save('RESTRICTEDMFNCOOPERATIVETARIFFPOLs.mat','RM');
        save('RESTRICTEDCOOPERATIVETARIFFPOLs.mat','R');
    end

    cd(root);
end

function polishCoop(r,v)
% Re-run a clean U for the three BAS folders using their own UM as guess.
    base = fullfile(r,'Results','Trade talks',v);
    for k = {'Nash_Bas','Fact','Free'}
        cd(fullfile(base,k{1}));
        load DATA.mat
        mycalculations
        S =  load('UNRESTRICTEDMFNCOOPERATIVETARIFFBASs.mat','UM');
        U =  mycooperativetariff(LAMBDABAS,-5,10,0,UM);
        save('UNRESTRICTEDCOOPERATIVETARIFFBASs.mat','U');
    end
    cd(r);
end

function cleanMFiles(f)
% Remove copied code to keep results dirs tidy.
    d = dir(fullfile(f,'**','*.m'));
    for i = 1:numel(d), delete(fullfile(d(i).folder,d(i).name)); end
end

function recalcRestrictedCoop(r,vs),src=fullfile(r,'Results','Trade talks',vs{2},'Free','RESTRICTEDCOOPERATIVETARIFFPOLs.mat');tgt=fullfile(r,'Results','Trade talks',vs{1},'Free');copyfile(src,tgt);cd(tgt);load DATA.mat;mycalculations;load RESTRICTEDCOOPERATIVETARIFFPOLs.mat;R=mymcooperativetariff(LAMBDAPOL,0,10,0,RESTRICTEDCOOPERATIVETARIFFPOLs);save RESTRICTEDCOOPERATIVETARIFFPOLs.mat R;cd(r);end
