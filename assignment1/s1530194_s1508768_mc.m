function [ aopt, fopt ] = s1530194_s1508768_mc(eval_budget)
% Calculates optimal power grid lay-out
%   Detailed explanation goes here
    bounds = load('para119.mat');    
    n=15;

    
    flag = 0;
    i=0;
    aopt = random_bounds(n,bounds);
    while flag==0
        i = i + 1;
        fprintf('In loop %d \n', i)
        aopt = random_bounds(n,bounds);
        flag = valid_119(aopt);
    end
    best_calc = calculation_119(aopt);
    fopt(1) = best_calc;
    reverseStr = '';
    for i = 2:eval_budget
        percentdone = 100*i/eval_budget;
        msg = sprintf('Percent done: %3.1f', percentdone);
        fprintf([reverseStr,msg]);
        reverseStr = repmat(sprintf('\b'),1,length(msg));
        
        a = random_bounds(n,bounds);
        flag = 0;
        while flag==0
            a = random_bounds(n,bounds);
            flag = valid_119(a);
        end

        if (calculation_119(a) < best_calc)
            aopt = a;
            best_calc = calculation_119(aopt);
        end
        fopt(i) = best_calc;
    end

end

function [a] = random_bounds(n,bounds)
    %Returns a random vector of length n which complies with the limits in
    %para119.mat

    for j = 1:n
        a(j) = randi([bounds.para.lb(j), bounds.para.ub(j)]);
    end
    a = a.';
end
        
        
        
        
        
        
        
        
        
        
