function [xopt, fopt] = s1508768_s1530194_es(fitnessfct, N, lb, ub, eval_budget)

mu = 2;
lambda = 10;
tau_prime = 1/sqrt(2*N);
tau = 1/sqrt(2*sqrt(N));

pa = rand(mu,N).*(repmat(ub,mu,1)-repmat(lb,mu,1)) + repmat(lb,mu,1);
si = 0.4*ones(mu,N);

eval_count = 0;
fopt_history = [];
while eval_count < eval_budget
    %decide which parents to pick (parents can be chosen multiple times)
    pick_parents = randi([1,mu],[1,lambda]);
    children = pa(pick_parents,:);
    si_children = si(pick_parents,:);

    si_children_mut = si_children .* exp(tau_prime * randn() + tau * randn(lambda,N));
    children_mut = children + si_children_mut .* randn(lambda,N);
    
    %select best individuals here using (mu,lambda) strategy
    for i = 1:lambda
        fitness(i) = fitnessfct(children_mut(i,:));
    end
    
    %select individuals
    [pa,idx] = select_tournament(children_mut, fitness);
    pa = pa(1:mu,:);
    %select sigmas
    si = si_children(idx,:);
    
    
    fopt_history = [fopt_history, fitness(1)];
    
    semilogy(fopt_history)
    xlim([1 eval_budget/lambda])     
    drawnow()
    
    eval_count = eval_count + lambda;
end




%{
nq := n(n-1)/2;
for i:=1 to n do
?u[i] := ?[i] * Ni(0,1);
for k:=1 to n-1 do
n1 := n-k;
n2 := n;
for i:=1 to k do
d1 := ?u[n1]; d2:= ?u[n2];
?c[n2] := d1*sin(?[nq])+ d2*cos(?[nq]);
?c[n1] := d1*cos(?[nq])- d2*sin(?[nq]);
n2 := n2-1;
nq := nq-1;
od
od
%}

end

function [a, idx] = select_tournament(P, f)
    %sort results by fitness and returns sorted population and the indicess that sort it
    [v,idx] = sort(f, 'ascend');
    a = P(idx,:);
end
