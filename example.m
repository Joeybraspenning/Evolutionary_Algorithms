
flag = 0;
tic;
while flag==0
    individual = randi([1,20],[15,1]);
    flag = valid_119(individual);
end
toc;

p_loss = calculation_119(individual);

