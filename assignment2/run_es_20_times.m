function [xopt_list, fopt_list] =  run_es_20_times(bbf)
    for loop_counter = 1:20
        [xopt, fopt] = s1508768_s1530194_es(bbf, 30, -100*ones(1,30), 100*ones(1,30), 10000);
        xopt_list(loop_counter,:) = xopt;
        fopt_list(loop_counter,:) = fopt; %save the fopt history
    end
    
    csvwrite(strcat(func2str(bbf), '_fopt_history.csv'), fopt_list);
    csvwrite(strcat(func2str(bbf), '_xopt.csv'), xopt_list);
end