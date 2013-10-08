% GA_SOLVER - Standard genetic algorithm solver
%
% Usage: ga_solver (p, n_gen)
%
% Where   p is a vector of motor performance parameters:
%         p = [sf eff pf Tb Tlr Ilr]
%           sf = full-load slip
%           eff = full-load efficiency
%           pf = full-load power factor
%           T_b = breakdown torque (as % of FL torque)
%           T_lr = locked rotor torque (as % of FL torque)
%           I_lr = locked rotor current
%		  n_gen is the maximum number of generations
%
% OUTPUT:   x is a vector of motor equivalent parameters:
%           x = [Rs Xs Xm Rr1 Xr1 Rr2 Xr2 Rc]

function [z err conv] = ga_solver(p, n_gen)

% Settings
pop = 20;     % population in each generation
n_r = 15;      % number of members retained for mating
n_e = 2;      % number of elite children per generation
c_f = 0.8;    % crossover fraction
% standard deviation weighting vector for mutation noise
sigma = [0.01 0.01 0.33 0.01 0.01 0.01 0.01 6.67]; 

w = [0.15 0.15 5 0.15 0.3 0.15 0.15 100]; % Weighting vector
n_gen = 30;    % maximum number of generations
err_tol = 0.00001; % error tolerance

% Human-readable motor performance parameters
% And base value initialisation
sf = p(1);                          % Full-load slip (pu)
eff = p(2);                         % Full-load efficiency (pu)
pf = p(3);                          % Full-load power factor (pu)
T_fl = pf * eff / (1 - sf);         % Full-load torque (pu)
T_b = p(4) * T_fl;                  % Breakdown torque (pu)
T_lr = p(5) * T_fl;                 % Locked rotor torque (pu)
i_lr = p(6);                        % Locked rotor current (pu)
Pm_fl = pf * eff;                   % Mechanical power (at FL)
Q_fl = sin(acos(pf));               % Full-load reactive power (pu)

% Formulate solution
pqt = [Pm_fl Q_fl T_b T_lr i_lr eff];

% Create initial population
x = (w' * ones(1,pop))' .* rand(pop,8);
gen = 1;
for i=1:pop    
    y = (pqt - calc_pqt(sf,x(i,:)))./pqt;
    err(i) = y*y';
    
    if err(i) < err_tol
        z = x(i,:);
        conv = 1;
        err = err(i);
        return
    end
end

% Genetic algorithm
for gen=2:n_gen
    
    % Select for fitness
    [fitness index] = sort(err',1);

    % Create next generation
    x_mate = x(index(1:n_r),:);         % select mating pool

    % Elite children (select best "n_e" children for next generation)
    x_new = x(index(1:n_e),:);

    % Crossover (random weighted average of parents)
    n_c = round((pop - n_e)*c_f);       % number of crossover children

    for j=1:n_c
        i_pair = ceil(n_r*rand(2,1));       % generate random pair of parents
        weight = rand(1,8);                    % generate random weighting
        % Crossover parents by weighted blend to generate new child
        x_new((n_e + j),:) = weight .* x_mate(i_pair(1),:) + (ones(1,8) - weight) .* 
			x_mate(i_pair(2),:);
    end

    % Mutation (gaussian noise added to parents)
    n_m = pop - n_e - n_c;       % number of mutation children

    for k=1:n_m
        % Select random parent from mating pool and add white noise
        x_new((n_e + n_c + k),:) = abs(x_mate(ceil(n_r*rand()),:) + sigma.*randn(1,8));
    end
    
    x = x_new;

    for i=1:pop    
        y = (pqt - calc_pqt(sf,x(i,:)))./pqt;
        err(i) = y*y';

        if err(i) < err_tol
            z = x(i,:);
            conv = 1;
            err = err(i);
            return
        end
    end
    
    % If the last generation, then output best results
    if gen == n_gen
        [fitness index] = sort(err',1);
        z = x(index(1),:);
        conv = 0;
        err = fitness(1);
    end
    
end

end
