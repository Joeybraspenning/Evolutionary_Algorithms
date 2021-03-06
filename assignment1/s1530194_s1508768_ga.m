function [opt, fopt] = s1530194_s1508768_ga(eval_budget)
% Calculates optimal power grid lay-out
%{
This function applies a genetic (7+28) algorithm to find the power grid lay out with the lowest power loss
The population is first initalized at random
As long as the evaluation budget is not reached, the population is crossed-over, then mutated and at last selected.
Some statistics are done to ensure that the best individual is captured in each generation
Input:
eval_budget = number of times the function calculate119 may be called
Output:
opt = optimal grid configuration for fopt
fopt = lowest powerloss / best fitness
The relevant parameters are:
n = vector length, number of loops in the power grid
mu = parent population size
offspring_ratio = number of children for each parent. Children population is of size mu*offspring_ratio
pc = cross-over probability
pm = mutation probability
%}
%Define parameters
  bounds = load('para119.mat');    
  n=15;
  mu = 7;
  offspring_ratio = 4;
  pc = 1/mu;
  pm = 2/n;    


  % Initialize population
  [F, P, f, evalcount, OP] = initialize(mu, n, bounds);

      % Evolution loop
  while evalcount < eval_budget
  [Pnew] = cross_over(pc, n, mu, P, offspring_ratio);
  [Pnew_new] = mutate(pm, n, Pnew, bounds, mu, offspring_ratio);
  [P, f] = selection(Pnew_new, P, mu, offspring_ratio, f);
  [F, evalcount, OP] = statistics(P, f, F, offspring_ratio, mu, evalcount, n, OP);
    
  end
  %Find the best individual of all generations and return the fitness and grid lay out
  [F_min, F_idx] = min(F);
  fopt = F_min;
  opt = OP(:, F_idx);
end


function [F, P, f, evalcount, OP] = initialize(mu, n, bounds)
    %{
    Input:
    global parameters mu, n, bounds
    evalcount = number of calls to calculate119 so far
    F = list of best fitness values for each generation
    Output:
    P = parent population
    f = fitness of parent population
    OP = list with the grid layout for the best fitness value for each generation
    
    Intialize the population with mu valid individuals and compute their fitness
    %}
    evalcount = 0;
    F = [];
    OP=[];
    for i = 1:mu
        %Initialize a random individual until a
        %valid solution is achieved
        j = 0;
        flag = 0;
        while flag==0
            j = j + 1;
            aopt = random_bounds(n,bounds);
            flag = valid_119(aopt);
        end
        P(:,i) = aopt;
        f(i) = calculation_119(P(:,i));
        evalcount = evalcount + 1;
    end
    [f_min, f_idx] = min(f);
    F = [F; f_min];
    OP = [OP, P(:,f_idx)];
end

function [Pnew] = cross_over(pc, n, mu, p, offspring_ratio)
    %{
    Perform cross-over between two randomly selected individuals from the population
    Cross-over is only performed with a probability of pc 
    
    Input:
    global parameters pc, n, mu, offspringratio
    p = current population before cross ver
    Output:
    Pnew = population after cross over
    %}
  for i = 1:offspring_ratio*mu
      if (rand() < pc)
        %select individual from population and attempt cross-over until a
        %valid solution is achieved
        flag = 0;
        j = 0;
        while flag==0
            j = j + 1;
            %select two random individuals
            index = randi([1,mu], [1,2]);
            %select the cross-over location at random
            crossloc = randi(n-1,1,1);
            
            Pnew(:,i) = [p(1:crossloc, index(1));p(crossloc:n-1,index(2))]; % crossover
            
            flag = valid_119(Pnew(:,i));
            if j > 1000
                Pnew(:,i) = p(:,index(1));
                break
            end
            %}
        end
      else
        %If cross-over does not occur the relevant individual is copied
        Pnew(:,i) = p(:,ceil(i/offspring_ratio)); % copy
      end
  end
end

function [Pnew_new] = mutate(pm, n, Pnew, bounds, mu, offspring_ratio)
  %{
  Perform mutation on random location of the individual with probability pm
  If a location is selected for mutation, a random number in bounds for
   that specific location will be chosen as the new value
  Input:
  global parameters pm, n, mu, bounds, offspring_ratio
  Pnew = population after cross over
  Output:
  Pnew_new = population after mutation
  %}
  for i = 1:offspring_ratio*mu
      %mutate individual from population until a valid solution is achieved
      flag = 0;
      while flag == 0
        %select random mutation location(s)
        mutloc = (rand(n,1) < pm);
        if sum(mutloc)>0
            Pnew_new(:,i) = mutate_in_bounds(Pnew(:,i), mutloc, bounds);% mutation
        else
            Pnew_new(:,i) = Pnew(:,i); % copy
        end
        flag = valid_119(Pnew_new(:,i));
      end
    end
end

function [P, f] = selection(Pnew_new, P, mu, offspring_ratio, f)
    %{
    Selects the fittest mu individuals from the combined parent + offspring population
    The selection method applied is tournament selection
    Input:
    global parameters u, offspring_ratio
    P = parent population
    Pnew_new = offspring after both cross-over and mutation
    f = fitness of parent population
    Output:
    P = new parent population
    f = fitness of new parent population
    %}
    for i = 1:offspring_ratio*mu
      f_new(i) = calculation_119(Pnew_new(:,i));
    end    
    
    % Combine the parent and offspring generations
    f_comb = [f,f_new];
    P_comb = [P,Pnew_new];
    %Sort the parent and offspring populations 
    [p_comb, idx] = select_tournament(P_comb, f_comb);
    
    %Pick the best mu individuals from the combined population
    P = p_comb(:,1:mu);
    f = f_comb(idx(1:mu));
end

function [F, evalcount, OP] = statistics(P, f, F, offspring_ratio, mu, evalcount, n, OP)
    %{
    Keeps track of the best individual and grid lay out for each generation
    Updates the evalcount to reflect the newest generation
    Input:
    global parameters offspring_ratio, mu, n
    P = new parent population after the function selection()
    f = fitness of new parent population
    F = list of best fitness values for each geneartion
    OP = list of best grid configurations for each generation
    evalcount = number of times the function calculate 119 is called so far
    Output:
    F = list of best fitness values for each geneartion, now appended with the value for the current generation
    evalcount = number of times the function calculate 119 is called so far, updates for the current generation
    OP = list of best grid configurations for each generation, appended with the current best configuration
    %}
    %Find best fitness for this generation
    [fopt, optindex] = min(f);
    F = [F; fopt];
    %Find best grid lay out for this generation
    opt = P(:,optindex);
    OP = [OP, opt];
    %Update the evalcount
    for i = 1:offspring_ratio*mu
      evalcount = evalcount + 1;
    end
end


function [a] = random_bounds(n,bounds)
    %Returns a random vector of length n which complies with the limits in
    %para119.mat
    %Input:
    %global parameters n, bounds
    %Output:
    %a = random vector of length n within bounds

    for j = 1:n
        a(j) = randi([bounds.para.lb(j), bounds.para.ub(j)]);
    end
    a = a.';
end

function [a, idx] = select_tournament(P, f)
    %{
    sort results by fitness and returns sorted population and the indicess that sort it
    Input:
    P = population to be selected
    f = fitness values of P
    Output:
    a = population P sorted by f
    idx = indicess that sort P
    %}
    [v,idx] = sort(f, 'ascend');
    a = P(:,idx);
end

function a = mutate_in_bounds(P, mutloc, bounds)
  %{
  Performs the actual mutation
  If a location in the individual is selected for mutation,
  that location is assigned a new random number drawn at random from within the bounds
  specified in the para119.mat file for this location
  The mutated vector is returned
  Input:
  global parameter bounds
  P = individual to be mutated
  mutloc = mutation locations
  Output:
  a = mutated individual
  %}
    for j =  1:length(mutloc)
        if mutloc(j)==1
            a(j) = randi([bounds.para.lb(j), bounds.para.ub(j)]);
        else 
            a(j) = P(j);
        end
    end
    a = a.';
end


