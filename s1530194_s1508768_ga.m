function [ opt, fopt ] = s1530194_s1508768_ga(eval_budget)
% Calculates optimal power grid lay-out
%   Detailed explanation goes here


%These settings generate a minimum powerloss of 787.0179
%    n=15;
% 	mu = 10;
%    pc = 1/mu;
%    pm = 2/n;    
%    evalcount = 0;
% with switch positions:[21,16,2,18,1,8,10,10,9,1,2,18,12,1,17]


    bounds = load('para119.mat');    
    n=15;
 	mu = 20;
    pc = 1/mu;
    pm = 6/n;    
    evalcount = 0;
    
    
    % Initialize population
    for i = 1:mu
        aopt = random_bounds(n,bounds);
        j = 0;
        flag = 0;
        while flag==0
            j = j + 1;
            aopt = random_bounds(n,bounds);
            flag = valid_119(aopt);
        end
        P(:,i) = aopt;
        f(i) = calculation_119(P(:,i));
    end
    
      % Evolution loop
  while evalcount < eval_budget
    [p, idx] = select_tournament(P, f);
    % Generate new population (recombination, mutation)
    for i = 1:mu
      %fprintf('In while %d forloop %d \n', evalcount, i)
      
      %adjust cross-over when a local minimum is reached
      
      %if evalcount == 2000
      %    pc = pc*2;
      %    pm = pm*1.5;
      %end
      %if evalcount > 2000
      %    p(:, floor(mu/4):floor(mu/2)) = p(:, floor(3*mu/4):mu);
      %end
      
      
      if (rand() < pc)
        %select individual from population and attempt cross-over until a
        %valid solution is achieved
        flag = 0;
        j = 0;
        while flag==0
            j = j + 1;
            crossloc = randi(n-1,1,1);
            index = randi([1,floor(mu/2)], [1,2]);
            Pnew(:,i) = [p(1:crossloc, index(1));p(crossloc:n-1,index(2))]; % crossover
            flag = valid_119(Pnew(:,i));
            if j > 1000
                %disp('No proper crossover found')
                Pnew(:,i) = p(:,index(1));
                break
            end
        end
      else
        Pnew(:,i) = p(:,i); % copy
      end
      
      %attempt mutation until a valid solution is found
      flag = 0;
      while flag == 0
        mutloc = (rand(n,1) < pm);
        if sum(mutloc)>0
            Pnew_new(:,i) = mutate_in_bounds(Pnew(:,i), mutloc, bounds);% mutation
        else
            Pnew_new(:,i) = Pnew(:,i); % copy
        end
        flag = valid_119(Pnew_new(:,i));
      end
    end

    % Decode and evaluate
    for i = 1:mu
      %g(:,i) = feval(decodefct, P(:,i));
      f_new(i) = calculation_119(Pnew_new(:,i));
    end    
    
    % Replace old population by new population, by choosing
    %the best individuals from both populations combined (comb)
    f_comb = [f,f_new];
    P_comb = [P,Pnew_new];
    [p_comb, idx] = select_tournament(P_comb, f_comb);
    P = p_comb(:,1:mu);
    f = f_comb(idx(1:mu));
    
    
    % Statistics administration
    [fopt, optindex] = min(f);
    opt = P(:,optindex);
    for i = 1:mu
      evalcount = evalcount + 1;
      histf(evalcount) = fopt;
    end
    
    fprintf('Best result: %10.4f kW\n', fopt)

    % Plot statistics
    clf
    subplot(2,1,1)
    semilogy(histf(1:evalcount))
    subplot(2,1,2)
    bar([1:n],opt)
    xlim([1 n])
    drawnow()
    
    %best result from literature: 869.7271 kW
  end
  %fprintf('Best result: %10.4f kW\n', fopt)
end


function [a] = random_bounds(n,bounds)
    %Returns a random vector of length n which complies with the limits in
    %para119.mat

    for j = 1:n
        a(j) = randi([bounds.para.lb(j), bounds.para.ub(j)]);
    end
    a = a.';
end

function [a, idx] = select_tournament(P, f)

    [v,idx] = sort(f, 'ascend');
    a = P(:,idx);
   
end

function a = mutate_in_bounds(P, mutloc, bounds)
    for j =  1:length(mutloc)
        if mutloc(j)==1
            a(j) = randi([bounds.para.lb(j), bounds.para.ub(j)]);
        else 
            a(j) = P(j);
        end
    end
    a = a.';
end

