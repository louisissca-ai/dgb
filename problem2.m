n = length(x);
opts = optimoptions('ga', 'Display', 'iter');

% 约束：sum(t) = 24
Aeq = ones(1, n);
beq = 24;

t_opt = ga(@(t) objFun(t, classicload, classicwind, classicsun, fenshijijia), ...
           n, [], [], Aeq, beq, zeros(n,1), ones(n,1), [], 1:n, opts);
t_opt = t_opt(:);
J_min = objFun(t_opt, classicload, classicwind, classicsun, fenshijijia);
fprintf('最小目标值 J = %.6g\n', J_min);
disp('最优 t 向量：'); disp(t_opt);

problem2_nh3load = t_opt*20.75*2;
problem2_load = classicload*6 + problem2_nh3load;
problem2_generate = classicwind*40 + classicsun*64;
problem2_buy = problem2_load - problem2_generate;
problem2_buy(problem2_buy<0) = 0;
problem2_sell(problem2_sell>0) = 0;
problem2_buycost = sum(problem2_buy .* fenshijijia);
problem2_dayload = sum(problem2_load);
problem2_daygenerate = sum(problem2_generate);
problem2_daybuy = sum(problem2_buy);
problem2_daysell = abs(sum(problem2_sell));
problem2_xnyzf = (problem1_dayload-problem1_daysell-problem1_daybuy)/problem1_daygenerate;
problem2_zydl = (problem1_daygenerate-problem1_daysell)/problem1_dayload;
problem2_xnysw = problem1_daysell/problem1_daygenerate;
fprintf('新能源自发自用电量占总可用发电量比例 = %.6g\n', problem2_xnyzf);
fprintf('总用电量中绿电占比 = %.6g\n', problem2_zydl);
fprintf('新能源上网电量比例 = %.6g\n', problem2_xnysw);

function J = objFun(t, classicload, classicwind, classicsun, fenshijijia)
    t = t(:);
    problem2_nh3load = t*20.75*2;
    problem2_load = classicload*6 + problem2_nh3load;
    problem2_generate = classicwind*40 + classicsun*64;
    problem2_buy = problem2_load - problem2_generate;
    problem2_sell = problem2_buy;
    problem2_buy(problem2_buy<0) = 0;
    problem2_sell(problem2_sell>0) = 0;
    problem2_buycost = sum(problem2_buy .* fenshijijia);
    problem2_sellcost = sum(problem2_sell*377.9);
    problem2_cost =   problem2_sellcost + problem2_buycost ...
                  + 0.75*2*sum(t)*2 ...
                  + 1.5*2*1000*0.2*60000/365/30 ...
                  + 10*2*sum(t)*100 ...
                  + 10*2*sum(t)*150 ...
                  + sum(classicsun)*64*1000*0.12 ...
                  + sum(classicwind)*40*1000*0.15;
    J = problem2_cost / (sum(t)*3);
end
