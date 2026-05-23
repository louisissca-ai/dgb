% problem4_2.m
% 问题四（2）：储能容量最优配置与有储能的最优调度
% 依赖 workspace 变量：wind1~wind6, sun1~sun4, classicload, fenshijijia

wind_vectors = {wind1(:), wind2(:), wind3(:), wind4(:), wind5(:), wind6(:)};
wind_names   = {'wind1', 'wind2', 'wind3', 'wind4', 'wind5', 'wind6'};
sun_vectors  = {sun1(:), sun2(:), sun3(:), sun4(:)};
sun_names    = {'sun1', 'sun2', 'sun3', 'sun4'};

classicload  = classicload(:);
fenshijijia  = fenshijijia(:);
n    = numel(wind_vectors{1});
lb_t = 0.1;

% =========================================================================
% 第一步：计算各场景弃电/缺电向量，确定最大弃电场景
% =========================================================================
fprintf('==================== 弃电/缺电向量分析 ====================\n');

slack_vectors  = cell(numel(wind_vectors), numel(sun_vectors));
best_max_slack = -Inf;
best_wind_idx  = 1;
best_sun_idx   = 1;

for wind_idx = 1:numel(wind_vectors)
    wind = wind_vectors{wind_idx};
    for sun_idx = 1:numel(sun_vectors)
        sun = sun_vectors{sun_idx};

        raw_ub = (wind*40 + sun*64 - classicload*6) / (20.75*2);
        ub_val = min(ones(n,1), raw_ub);

        slack = zeros(n,1);
        for i = 1:n
            if ub_val(i) < lb_t
                slack(i) = ub_val(i) - lb_t;   % 缺电：负值
            elseif raw_ub(i) > 1
                slack(i) = raw_ub(i) - 1;       % 弃电：正值
            end
        end

        slack_vectors{wind_idx, sun_idx} = slack;
        fprintf('  %s + %s 弃电/缺电向量：\n', wind_names{wind_idx}, sun_names{sun_idx});
        fprintf('    ['); fprintf(' %8.4g', slack); fprintf(' ]\n');

        if max(slack) > best_max_slack
            best_max_slack = max(slack);
            best_wind_idx  = wind_idx;
            best_sun_idx   = sun_idx;
        end
    end
end

fprintf('最大弃电场景：%s + %s（最大弃电量 = %.4g）\n', ...
        wind_names{best_wind_idx}, sun_names{best_sun_idx}, best_max_slack);

% =========================================================================
% 第二步：由最大弃电场景计算储能容量
% =========================================================================
best_slack = slack_vectors{best_wind_idx, best_sun_idx};
b = cumsum(best_slack);

[max_b, max_b_idx] = max(b);
min_b_after = min(b(max_b_idx:end));   % max(b) 之后的最小值
E = (max_b - min_b_after) * 20.75 * 2;

fprintf('b 向量：max = %.4g（位置 %d），max 之后最小 = %.4g\n', max_b, max_b_idx, min_b_after);
fprintf('储能容量 E = (%.4g - %.4g) * 20.75 * 2 = %.4g MW·h\n', ...
        max_b, min_b_after, E);

figure('Name', 'problem4-2：最大弃电场景储能分析', 'Color', 'w');
subplot(2,1,1);
bar(best_slack, 'FaceColor', [0.2 0.5 0.8]);
grid on;
title(sprintf('%s + %s：弃电(+)/缺电(-) 向量', ...
      wind_names{best_wind_idx}, sun_names{best_sun_idx}), 'Interpreter', 'none');
xlabel('时段'); ylabel('弃电(+) / 缺电(-)');

subplot(2,1,2);
plot(b, 'b-o', 'LineWidth', 1.2); grid on;
yline(max_b,       'r--', sprintf('max = %.4g（位置 %d）', max_b, max_b_idx));
yline(min_b_after, 'g--', sprintf('max 后最小 = %.4g', min_b_after));
title(sprintf('累积向量 b，储能容量 E = %.4g MW·h', E), 'Interpreter', 'none');
xlabel('时段'); ylabel('累积弃电量');

% =========================================================================
% 第三步：有储能的最优调度（24 场景）
% =========================================================================
fprintf('\n\n==================== 有储能最优调度（24 场景）====================\n');
fprintf('储能容量 E = %.4g MW·h\n\n', E);

opts_lp = optimoptions('linprog', 'Display', 'off');
storage_results = struct([]);
res_idx = 0;

for wind_idx = 1:numel(wind_vectors)
    wind = wind_vectors{wind_idx};
    wind_name = wind_names{wind_idx};

    for sun_idx = 1:numel(sun_vectors)
        sun = sun_vectors{sun_idx};
        sun_name = sun_names{sun_idx};

        fprintf('\n==================== 情况 %02d：%s + %s ====================\n', ...
                (wind_idx-1)*numel(sun_vectors)+sun_idx, wind_name, sun_name);

        % 决策变量 x = [t(n); SOC(n); ch(n); dis(n); buy(n); sell(n)]
        % 目标：最大化 sum(t) ↔ 最小化 -sum(t)
        f_lp  = [-ones(n,1); zeros(5*n,1)];
        lb_lp = [lb_t*ones(n,1); zeros(5*n,1)];
        ub_lp = [ones(n,1); E*ones(n,1); Inf(4*n,1)];

        % 等式约束 1：逐时段功率平衡
        % 20.75*2*t[i] + ch[i] - dis[i] - buy[i] + sell[i]
        %   = wind[i]*40 + sun[i]*64 - classicload[i]*6
        Aeq1 = [diag(repmat(20.75*2, n, 1)), zeros(n,n), ...
                eye(n), -eye(n), -eye(n), eye(n)];
        beq1 = wind*40 + sun*64 - classicload*6;

        % 等式约束 2：SOC 动态（初始 SOC[0] = 0）
        % SOC[i] - SOC[i-1] - ch[i] + dis[i] = 0
        SOC_coef = eye(n) - diag(ones(n-1,1), -1);
        Aeq2 = [zeros(n,n), SOC_coef, -eye(n), eye(n), zeros(n,n), zeros(n,n)];
        beq2 = zeros(n,1);

        Aeq = [Aeq1; Aeq2];
        beq = [beq1; beq2];

        [x_sol, ~, exitflag] = linprog(f_lp, [], [], Aeq, beq, lb_lp, ub_lp, opts_lp);

        if exitflag ~= 1
            fprintf('LP 求解失败（exitflag=%d），跳过。\n', exitflag);
            continue;
        end

        t_opt    = x_sol(1:n);
        soc_opt  = x_sol(n+1:2*n);
        ch_opt   = x_sol(2*n+1:3*n);
        dis_opt  = x_sol(3*n+1:4*n);
        buy_opt  = x_sol(4*n+1:5*n);
        sell_opt = x_sol(5*n+1:6*n);

        sum_t = sum(t_opt);
        fprintf('有储能 sum(t) = %.6g\n', sum_t);

        % 能量指标
        nh3load_vec = t_opt * 20.75 * 2;
        load_vec    = classicload*6 + nh3load_vec;
        gen_vec     = wind*40 + sun*64;

        dayload     = sum(load_vec);
        daygenerate = sum(gen_vec);
        daybuy      = sum(buy_opt);
        daysell     = sum(sell_opt);

        xnyzf = (dayload - daysell - daybuy) / daygenerate;
        zydl  = (daygenerate - daysell) / dayload;
        xnysw = daysell / daygenerate;

        % 成本计算（与 objFun 保持一致）
        buycost  = sum(buy_opt .* fenshijijia);
        sellcost = -sum(sell_opt) * 377.9;         % 卖电收益（负成本）
        cost_total = buycost + sellcost ...
                   + 0.75*2*sum_t*2 ...
                   + 1.5*2*1000*0.2*60000/365/30 ...
                   + 10*2*sum_t*100 ...
                   + 10*2*sum_t*150 ...
                   + sum(sun)*64*1000*0.12 ...
                   + sum(wind)*40*1000*0.15;

        J_storage       = cost_total / (sum_t * 3);
        annual_nh3      = sum_t * 3 * 365;          % 全年制氨量（吨）
        cap_utilization = sum_t / n;                % 年平均产能利用率

        fprintf('吨氨成本 J = %.6g\n',    J_storage);
        fprintf('全年制氨总量 = %.6g 吨\n', annual_nh3);
        fprintf('年平均产能利用率 = %.4f\n', cap_utilization);
        fprintf('xnyzf=%.6g，zydl=%.6g，xnysw=%.6g\n', xnyzf, zydl, xnysw);

        res_idx = res_idx + 1;
        storage_results(res_idx).case_idx        = res_idx;
        storage_results(res_idx).wind_idx        = wind_idx;
        storage_results(res_idx).wind_name       = wind_name;
        storage_results(res_idx).sun_idx         = sun_idx;
        storage_results(res_idx).sun_name        = sun_name;
        storage_results(res_idx).t_opt           = t_opt;
        storage_results(res_idx).soc_opt         = soc_opt;
        storage_results(res_idx).ch_opt          = ch_opt;
        storage_results(res_idx).dis_opt         = dis_opt;
        storage_results(res_idx).buy             = buy_opt;
        storage_results(res_idx).sell            = sell_opt;
        storage_results(res_idx).sum_t           = sum_t;
        storage_results(res_idx).daybuy          = daybuy;
        storage_results(res_idx).daysell         = daysell;
        storage_results(res_idx).dayload         = dayload;
        storage_results(res_idx).daygenerate     = daygenerate;
        storage_results(res_idx).xnyzf           = xnyzf;
        storage_results(res_idx).zydl            = zydl;
        storage_results(res_idx).xnysw           = xnysw;
        storage_results(res_idx).J               = J_storage;
        storage_results(res_idx).annual_nh3      = annual_nh3;
        storage_results(res_idx).cap_utilization = cap_utilization;
    end
end

% =========================================================================
% 第四步：绘图
% =========================================================================
plotStorageResults(storage_results, E);

% =========================================================================
% 局部函数
% =========================================================================
function plotStorageResults(storage_results, E)
    if isempty(storage_results)
        warning('没有可绘制的结果。');
        return;
    end

    sort_keys = [[storage_results.wind_idx]', [storage_results.sun_idx]'];
    [~, order] = sortrows(sort_keys, [1 2]);
    storage_results = storage_results(order);

    combo_labels = arrayfun(@(r) sprintf('%s+%s', r.wind_name, r.sun_name), ...
                            storage_results, 'UniformOutput', false);

    metric_fields = {'sum_t', 'J', 'annual_nh3', 'cap_utilization', 'xnyzf', 'zydl'};
    metric_names  = {'sum(t)', '吨氨成本(元/吨)', '全年制氨量(吨)', ...
                     '产能利用率', 'xnyzf', 'zydl'};

    figure('Name', sprintf('problem4-2：有储能(E=%.4gMWh)各场景指标', E), 'Color', 'w');

    for mi = 1:numel(metric_fields)
        values = arrayfun(@(r) r.(metric_fields{mi}), storage_results);
        valid  = values(~isnan(values));

        subplot(2, 3, mi);
        bar(values, 'FaceColor', [0.4 0.2 0.7]);
        grid on; hold on;

        mean_v = mean(valid);
        plot(xlim, [mean_v mean_v], 'r--', 'LineWidth', 1.2);

        [~, min_idx] = min(values);
        [~, max_idx] = max(values);
        scatter(min_idx, min(valid), 45, 'g', 'filled');
        scatter(max_idx, max(valid), 45, 'r', 'filled');

        title({metric_names{mi}, sprintf('均值=%.4g，最小=%.4g，最大=%.4g', ...
               mean_v, min(valid), max(valid))}, 'Interpreter', 'none');
        xlabel('wind+sun 组合');
        ylabel(metric_names{mi}, 'Interpreter', 'none');
        set(gca, 'XTick', 1:numel(combo_labels), ...
                 'XTickLabel', combo_labels, 'XTickLabelRotation', 45);
        xlim([0.5, numel(combo_labels)+0.5]);
        legend({'组合值','均值','最小','最大'}, 'Location','best','FontSize',8);
        hold off;
    end
end
