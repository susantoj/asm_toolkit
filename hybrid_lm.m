% HYBRID_LM - Hybrid genetic algorithm and LM solver to optimise Rs and Xr2
%             Uses LM_SOLVER1a as the base newton-raphson solver
%			  Note that LM_SOLVER1a is simply the conventional LM_SOLVER 
%			  (error term adjustment) with Rs and Xr2 as known inputs.
%             Optimises Rs and Xr2 with genetic algorithm
%

function [z err conv] = hybrid_lm(p, max_iter)

% Settings
pop = 15;     % population in each generation
n_r = 10;      % number of members retained for mating
n_e = 2;      % number of elite children per generation
c_f = 0.8;    % crossover fraction
sigma = 0.01; % standard deviation for mutation noise
n_gen = 10;    % maximum number of generations

% Create initial population
Rs = 0.15.*rand(pop,2)
gen = 1
for i=1:pop    
    [x(i,:) iter(i) err(i) conv(i)] = lm_solver1a(p, Rs(i,1), Rs(i,2), 1e-7, 5, max_iter);
    if conv(i) == 1
        z = x(i,:);
        conv = 1;
        err = err(i);
        return
    end
end

for gen=2:n_gen
    
    gen
    
    % Select for fitness
    [fitness index] = sort(err',1);

    % Create next generation
    Rs_mate = Rs(index(1:n_r),:);         % select mating pool

    % Elite children (select best "n_e" children for next generation)
    Rs_new = Rs(index(1:n_e),:);

    % Crossover (random weighted average of parents)
    n_c = round((pop - n_e)*c_f);       % number of crossover children

    for j=1:n_c
        i_pair = ceil(n_r*rand(2,1));       % generate random pair of parents
        weight = rand();                    % generate random weighting
        % Crossover parents by weighted blend to generate new child
        Rs_new((n_e + j),:) = weight .* Rs_mate(i_pair(1),:) + (1 - weight) .* 
			Rs_mate(i_pair(2),:);
    end

    % Mutation (gaussian noise added to parents)
    n_m = pop - n_e - n_c;       % number of mutation children

    for k=1:n_m
        % Select random parent from mating pool and add white noise
        Rs_new((n_e + n_c + k),:) = abs(Rs_mate(ceil(n_r*rand()),:) + sigma*randn(1,2));
    end
    
    Rs = Rs_new

    for i=1:pop
    [x(i,:) iter(i) err(i) conv(i)] = lm_solver1a(p, Rs(i,1), Rs(i,2), 1e-7, 5, max_iter);
        if conv(i) == 1
            z = x(i,:);
            conv = 1;
            err = err(i);
            return
        end
    end
    
    % If the last generation, then output best results
    if gen == n_gen
        [fitness index] = sort(err',1);
        z = x(index(1),:)
        conv = 0
        err = fitness(1)
        Rs_best = Rs(index(1),1)
        Xr2_best = Rs(index(1),2)
    end  
end

end
