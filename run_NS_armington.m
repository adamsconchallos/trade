par.N=3; par.S=2;
par.a = [1.0  1.3;
         1.1  1.1;
         0.9  1.4];        % N×S
par.L = [1; 1; 1];
par.sigma = [4; 6];        % S×1
par.mu = [0.6 0.4;
          0.6 0.4;
          0.6 0.4];        % rows sum to 1
par.kappa = ones(3,3,2);   % start frictionless
par.tau   = zeros(3,3,2);  % start free trade
sol = solve_GE_NS(par);     % solves for w (with w3=1 as numeraire)
disp(sol.out.w), disp(sol.out.W)
%%

par.sigma=4; par.aH=1; par.aF=1.2; par.LH=1; par.LF=1;
T = compare_tariff_cases(par, [0 1; 0 1]);
disp(T)
