function [xopt, fopt] = s1508768_s1530194_es(fitnessfct, N, lb, ub, eval_budget)
    %{
    fitnessfct = fitnessfunction to be optimised
    N = dimensionality of the problem
    lb = lower bounds for each dimension
    ub = upper bounds for each dimension
    eval_budget = maximum number of allowed evaluations of fitnessfct
    %}


    mu = 5;
    lambda = 25;
    tau_prime = 1/sqrt(2*N);
    tau = 1/sqrt(2*sqrt(N));
    beta = 0.0873; %rotation constant corresponds to 5 degrees

    pa = rand(mu,N).*(repmat(ub,mu,1)-repmat(lb,mu,1)) + repmat(lb,mu,1);
    si = 0.4*ones(mu,N); %was 0.4
    alpha = rand(mu,N) * 2 * pi - pi;

    eval_count = 0;
    fopt_history = [];
    xopt_history = [];
    while eval_count < eval_budget
        %decide which parents to pick (parents can be chosen multiple times)
        pick_parents = randi([1,mu],[1,lambda]);
        children = pa(pick_parents,:);
        si_children = si(pick_parents,:);
        alpha_children = alpha(pick_parents,:);

        si_children_mut = si_children .* exp(tau_prime * randn() + tau * randn(lambda,N));

        alpha_children_mut = alpha_children + beta * randn(lambda,N);

        %make sure that the angles do not come outside of the bounds -pi, piz
        if sum(abs(alpha_children_mut) > pi) > 0
            alpha_children_mut(abs(alpha_children_mut) > pi) = alpha_children_mut(alpha_children_mut > pi) - 2*pi * sign(alpha_children_mut(alpha_children_mut > pi));
        end

        delta_x = si_children_mut .* randn(lambda,N);
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

        children_mut = children + delta_x;
        %children_mut = children + si_children_mut .* randn(lambda,N);

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

        %select best individuals here using (mu,lambda) strategy
        for i = 1:lambda
            fitness(i) = fitnessfct(children_mut(i,:));
        end

        %select individuals
        [pa,idx] = select_tournament(children_mut, fitness);
        pa = pa(1:mu,:);
        %select sigmas
        si = si_children_mut(idx,:);
        alpha = alpha_children_mut(idx,:);


        fopt_history = [fopt_history, fitness(1)];
        xopt_history = [xopt_history; pa(1,:)];

        semilogy(fopt_history)
        xlim([1 eval_budget/lambda])     
        drawnow()

        eval_count = eval_count + lambda;
    end
    
    [value, minloc] = min(fopt_history);
    fopt = fopt_history;
    xopt = xopt_history(minloc,:);

end

function [a, idx] = select_tournament(P, f)
    %sort results by fitness and returns sorted population and the indicess that sort it
    [v,idx] = sort(f, 'ascend');
    a = P(idx,:);
end
