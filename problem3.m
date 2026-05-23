wind_vectors = {wind1(:), wind2(:), wind3(:), wind4(:), wind5(:), wind6(:)};
wind_names = {'wind1', 'wind2', 'wind3', 'wind4', 'wind5', 'wind6'};
sun_vectors = {sun1(:), sun2(:), sun3(:), sun4(:)};
sun_names = {'sun1', 'sun2', 'sun3', 'sun4'};

classicload = classicload(:);
fenshijijia = fenshijijia(:);
n = numel(wind_vectors{1});

% t 为连续变量，范围 [0.1, 1]，使用 fmincon (SQP) 求解
opts = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');

% 分别在以下约束下运行优化：sum(t) = beq
beq_values = [24, 21, 18, 15, 12];
Aeq = ones(1, n);
lb = 0.1 * ones(n, 1);
ub = ones(n, 1);

problem3_results = struct([]);
result_idx = 0;

for wind_idx = 1:numel(wind_vectors)
    wind = wind_vectors{wind_idx};
    wind_name = wind_names{wind_idx};

    for sun_idx = 1:numel(sun_vectors)
        sun = sun_vectors{sun_idx};
        sun_name = sun_names{sun_idx};

        fprintf('\n\n==================== 情况 %02d：%s + %s ====================\n', ...
                (wind_idx - 1)*numel(sun_vectors) + sun_idx, wind_name, sun_name);

        for beq_idx = 1:numel(beq_values)
            beq = beq_values(beq_idx);

            fprintf('\n========== %s + %s, beq = %g, sum(t) = %g ==========\n', ...
                    wind_name, sun_name, beq, beq);

            % 初始点：均匀分布并截断到 [0.1, 1]，满足等式约束
            x0 = (beq / n) * ones(n, 1);

            t_opt = fmincon(@(t) objFun(t, classicload, wind, sun, fenshijijia), ...
                            x0, [], [], Aeq, beq, lb, ub, [], opts);
            t_opt = t_opt(:);
            J_min = objFun(t_opt, classicload, wind, sun, fenshijijia);
            fprintf('最小目标值 J = %.6g\n', J_min);
            disp('最优 t 向量：'); disp(t_opt);

            metrics = calcProblem3Metrics(t_opt, classicload, wind, sun, fenshijijia);

            problem3_nh3load = metrics.nh3load;
            problem3_load = metrics.load;
            problem3_generate = metrics.generate;
            problem3_buy = metrics.buy;
            problem3_sell = metrics.sell;
            problem3_buycost = metrics.buycost;
            problem3_dayload = metrics.dayload;
            problem3_daygenerate = metrics.daygenerate;
            problem3_daybuy = metrics.daybuy;
            problem3_daysell = metrics.daysell;
            problem3_xnyzf = metrics.xnyzf;
            problem3_zydl = metrics.zydl;
            problem3_xnysw = metrics.xnysw;

            fprintf('新能源自发自用电量占总可用发电量比例 = %.6g\n', problem3_xnyzf);
            fprintf('总用电量中绿电占比 = %.6g\n', problem3_zydl);
            fprintf('新能源上网电量比例 = %.6g\n', problem3_xnysw);

            result_idx = result_idx + 1;
            problem3_results(result_idx).case_idx = result_idx;
            problem3_results(result_idx).wind_idx = wind_idx;
            problem3_results(result_idx).wind_name = wind_name;
            problem3_results(result_idx).sun_idx = sun_idx;
            problem3_results(result_idx).sun_name = sun_name;
            problem3_results(result_idx).beq = beq;
            problem3_results(result_idx).t_opt = t_opt;
            problem3_results(result_idx).J_min = J_min;
            problem3_results(result_idx).nh3load = problem3_nh3load;
            problem3_results(result_idx).load = problem3_load;
            problem3_results(result_idx).generate = problem3_generate;
            problem3_results(result_idx).buy = problem3_buy;
            problem3_results(result_idx).sell = problem3_sell;
            problem3_results(result_idx).buycost = problem3_buycost;
            problem3_results(result_idx).dayload = problem3_dayload;
            problem3_results(result_idx).daygenerate = problem3_daygenerate;
            problem3_results(result_idx).daybuy = problem3_daybuy;
            problem3_results(result_idx).daysell = problem3_daysell;
            problem3_results(result_idx).xnyzf = problem3_xnyzf;
            problem3_results(result_idx).zydl = problem3_zydl;
            problem3_results(result_idx).xnysw = problem3_xnysw;
        end
    end
end

plotProblem3Distributions(problem3_results, beq_values);

function metrics = calcProblem3Metrics(t, classicload, wind, sun, fenshijijia)
    t = t(:);
    classicload = classicload(:);
    wind = wind(:);
    sun = sun(:);
    fenshijijia = fenshijijia(:);

    metrics.nh3load = t*20.75*2;
    metrics.load = classicload*6 + metrics.nh3load;
    metrics.generate = wind*40 + sun*64;
    metrics.buy = metrics.load - metrics.generate;
    metrics.sell = metrics.buy;
    metrics.buy(metrics.buy<0) = 0;
    metrics.sell(metrics.sell>0) = 0;
    metrics.buycost = sum(metrics.buy .* fenshijijia);
    metrics.dayload = sum(metrics.load);
    metrics.daygenerate = sum(metrics.generate);
    metrics.daybuy = sum(metrics.buy);
    metrics.daysell = abs(sum(metrics.sell));

    metrics.xnyzf = (metrics.dayload-metrics.daysell-metrics.daybuy)/metrics.daygenerate;
    metrics.zydl = (metrics.daygenerate-metrics.daysell)/metrics.dayload;
    metrics.xnysw = metrics.daysell/metrics.daygenerate;
end

function J = objFun(t, classicload, wind, sun, fenshijijia)
    t = t(:);
    classicload = classicload(:);
    wind = wind(:);
    sun = sun(:);
    fenshijijia = fenshijijia(:);

    problem3_nh3load = t*20.75*2;
    problem3_load = classicload*6 + problem3_nh3load;
    problem3_generate = wind*40 + sun*64;
    problem3_buy = problem3_load - problem3_generate;
    problem3_sell = problem3_buy;
    problem3_buy(problem3_buy<0) = 0;
    problem3_sell(problem3_sell>0) = 0;
    problem3_buycost = sum(problem3_buy .* fenshijijia);
    problem3_sellcost = sum(problem3_sell*377.9);
    problem3_cost =   problem3_sellcost + problem3_buycost ...
                  + 0.75*2*sum(t)*2 ...
                  + 1.5*2*1000*0.2*60000/365/30 ...
                  + 10*2*sum(t)*100 ...
                  + 10*2*sum(t)*150 ...
                  + sum(sun)*64*1000*0.12 ...
                  + sum(wind)*40*1000*0.15;
    J = problem3_cost / (sum(t)*3);
end

function plotProblem3Distributions(problem3_results, beq_values)
    metric_fields = {'xnyzf', 'zydl', 'xnysw', 'J_min', 'daybuy', 'daysell'};
    metric_names = {'problem3\_xnyzf', 'problem3\_zydl', 'problem3\_xnysw', ...
                    'J', 'problem3\_daybuy', 'problem3\_daysell'};

    for beq_idx = 1:numel(beq_values)
        beq = beq_values(beq_idx);
        beq_mask = [problem3_results.beq] == beq;
        beq_results = problem3_results(beq_mask);

        if isempty(beq_results)
            warning('beq = %g 没有可绘制的结果。', beq);
            continue;
        end

        sort_keys = [[beq_results.wind_idx]', [beq_results.sun_idx]'];
        [~, order] = sortrows(sort_keys, [1 2]);
        beq_results = beq_results(order);

        combo_labels = arrayfun(@(r) sprintf('%s+%s', r.wind_name, r.sun_name), ...
                                beq_results, 'UniformOutput', false);

        figure('Name', sprintf('beq=%g 下 24 种 wind/sun 组合指标分布', beq), ...
               'Color', 'w');

        for metric_idx = 1:numel(metric_fields)
            metric_field = metric_fields{metric_idx};
            metric_name = metric_names{metric_idx};
            values = arrayfun(@(r) r.(metric_field), beq_results);

            subplot(2, 3, metric_idx);
            bar(values, 'FaceColor', [0.2 0.45 0.75]);
            grid on;
            hold on;

            valid_values = values(~isnan(values));
            mean_value = mean(valid_values);
            min_value = min(valid_values);
            max_value = max(valid_values);
            [~, min_idx] = min(values);
            [~, max_idx] = max(values);

            plot(xlim, [mean_value mean_value], 'r--', 'LineWidth', 1.2);
            scatter(min_idx, min_value, 45, 'g', 'filled');
            scatter(max_idx, max_value, 45, 'r', 'filled');

            title({sprintf('%s，beq = %g', metric_name, beq), ...
                   sprintf('均值=%.4g，最小=%.4g，最大=%.4g', ...
                           mean_value, min_value, max_value)}, ...
                  'Interpreter', 'none');
            xlabel('wind + sun 组合');
            ylabel(metric_name, 'Interpreter', 'none');
            set(gca, 'XTick', 1:numel(combo_labels), ...
                     'XTickLabel', combo_labels, ...
                     'XTickLabelRotation', 45);
            xlim([0.5, numel(combo_labels) + 0.5]);

            if min_value == max_value
                ylim([min_value - 0.5, max_value + 0.5]);
            end

            legend({'组合值', '均值', '最小值', '最大值'}, ...
                   'Location', 'best', 'FontSize', 8);
            hold off;
        end
    end
end
