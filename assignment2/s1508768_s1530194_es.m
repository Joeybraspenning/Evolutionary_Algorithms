function [xopt, fopt] = s1508768_s1530194_es(fitnessfct, N, lb, ub, eval_budget)
    %{
    fitnessfct = fitnessfunction to be optimised
    N = dimensionality of the problem
    lb = lower bounds for each dimension
    ub = upper bounds for each dimension
    eval_budget = maximum number of allowed evaluations of fitnessfct
    %}

    %Initialise global parameters
    mu = 3; %Parent population size
    lambda = 21; %Offspring size
    tau_prime = 1/sqrt(2*N); %constant for mutation
    tau = 1/sqrt(2*sqrt(N)); %constand for mutation
    beta = 0.0873; %rotation constant corresponds to 5 degrees
    
    %Intialise a random parent population
    pa = rand(mu,N).*(repmat(ub,mu,1)-repmat(lb,mu,1)) + repmat(lb,mu,1);
    si = 0.4*ones(mu,N);
    alpha = rand(mu,N) * 2 * pi - pi;
    
    %Calculate fitness of parent population
    for j = 1:mu
        fitness_old(j) = fitnessfct(pa(j,:));
    end
    %Book keeping
    eval_count = mu;
    
    %Intialise book keeping tools
    fopt_history = [];
    xopt_history = [];
    
    while eval_count < eval_budget
        %Create children
        pick_parents = randi([1,mu],[1,lambda]);
        children = pa(pick_parents,:);
        si_children = si(pick_parents,:);
        alpha_children = alpha(pick_parents,:);
        
        %Mutate children
        [children_mut, si_children_mut, alpha_children_mut] = mutate(children, si_children, alpha_children, tau, tau_prime, beta, lambda, N, lb, ub);

        %Calculate fitness for all children
        for i = 1:lambda
            fitness(i) = fitnessfct(children_mut(i,:));
        end

        %Select (mu+lambda)
        [pa, si, alpha, fitness_old] = selection(fitness, fitness_old, si_children_mut, si, alpha_children_mut, alpha, children_mut, pa, mu, lambda, N);

        %Do some book keeping
        fopt_history = [fopt_history, fitness_old(1)];
        xopt_history = [xopt_history; pa(1,:)];
        eval_count = eval_count + lambda;
    end
    
    %Return best fitness value and the best individual
    [value, minloc] = min(fopt_history);
    fopt = value;
    xopt = xopt_history(minloc,:);
end

function [a, idx] = select_tournament(P, f)
    %{
    Sort the population according to fitness value
    INPUT:
    P = population to be sorted
    f = fitness values of population P

    OUTPUT:
    a = sorted population
    idx = indexing that sorts populationb
    %}
    
    [v,idx] = sort(f, 'ascend');
    a = P(idx,:);
end

function [mut_individual, mut_si, mut_alpha] = mutate(children, si_children, alpha_children, tau, tau_prime, beta, lambda, N, lb, ub)
        %{
        Mutates the individuals, sigma values and rotation angles
        INPUT:
        Global parameters: tau, tau_prime, beta, lambda, N, lb, ub
        children = The children population, copies of parent population
        si_children = Sigma values for all dimensions of all children
        alpha_children = Rotation angles of children
        
        OUTPUT:
        mut_individual = The children after correlated mutation is applied
        mut_si = The mutated sigma values of the mut_individual
        mut_alpha = The mutated alpha of the mut_individual
        %}

        %mutate sigma values
        si_children_mut = si_children .* exp(tau_prime * randn() + tau * randn(lambda,N));
        
        %mutate rotation angle
        alpha_children_mut = alpha_children + beta * randn(lambda,N);

        %make sure that the angles do not come outside of the bounds -pi, piz
        if lambda>1
            if sum(abs(alpha_children_mut) > pi) > 0
                alpha_children_mut(abs(alpha_children_mut) > pi) = alpha_children_mut(alpha_children_mut > pi) - 2*pi * sign(alpha_children_mut(alpha_children_mut > pi));
            end
        else
            if abs(alpha_children_mut) > pi
                alpha_children_mut = alpha_children_mut - 2*pi * sign(alpha_children_mut)
            end
        end
        
        %The change to the children before correlated mutation
        delta_x = si_children_mut .* randn(lambda,N);
        
        %Apply correlated mutation = 2D rotation of the sigma values
        for m = 1:lambda
            for k = 1:N-1
                n1 = N-k;
                n2 = N;
                for j = 1:k
                    d1 = delta_x(m,n1);
                    d2 = delta_x(m,n2);
                    angle = alpha_children_mut(m,n2) - alpha_children_mut(m,n1);
                    delta_x(m,n2) = d1 * sin(angle) + d2 * cos(angle);
                    delta_x(m,n1) = d1 * cos(angle) - d2 * sin(angle);
                    n2 = n2-1;
                end
            end
        end 
        
        %Mutate the children
        children_mut = children + delta_x;

        %prevent the mutations from going beyond the bounds
        for i = 1:lambda
            for j = 1:N
                if children_mut(i,j) < lb(j)
                    children_mut(i,j) = lb(j);
                end
                if children_mut(i,j) > ub(j)
                    children_mut(i,j) = ub(j);
                end
            end
        end
        
        %Prepare the return values
        mut_individual = children_mut;
        mut_si = si_children_mut;
        mut_alpha = alpha_children_mut;
end

function [pa, si, alpha, fitness_old] = selection(fitness, fitness_old, si_children_mut, si, alpha_children_mut, alpha, children_mut, pa, mu, lambda, N)
        %{
        Selects the best mu individuals for (mu+lambda) selection
        INPUT: 
        Global parameters: mu, lambda, N
        fitness, fitness_old = fitness values of children and parents
        si_children_mut, si = sigmas of children and parents
        alpha_children_mut, alpha = rotation angles of children and parents
        children_mut, pa = children and parent individuals
        
        OUTPUT:
        pa = New parent population (best mu individuals)
        si = Sigmas of new parent population
        alpha = Rotation angles of new parent population
        fitness_old = fitness values of new parent population
        %}

        %Combine parent and child population
        fitness_all = [fitness, fitness_old];
        si_all = [si_children_mut; si];
        alpha_all = [alpha_children_mut; alpha];
        all_individuals = [children_mut; pa];
        
        %select individuals
        [pa,idx] = select_tournament(all_individuals, fitness_all);
        
        %Pick best mu individuals
        pa = pa(1:mu,:);
        
        %Pick sigmas of best mu individuals
        si_all = si_all(idx,:);
        si = si_all(1:mu,:);
        
        %Pick alphas of best mu individuals
        alpha_all = alpha_all(idx,:);
        alpha = alpha_all(1:mu, :);
        
        %Pick fitness of best mu individuals
        fitness_all = fitness_all(:, idx);
        fitness_old = fitness_all(:, 1:mu);
end
        
