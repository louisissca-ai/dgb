wind_vectors = {wind1(:), wind2(:), wind3(:), wind4(:), wind5(:), wind6(:)};
wind_names = {'wind1', 'wind2', 'wind3', 'wind4', 'wind5', 'wind6'};
sun_vectors = {sun1(:), sun2(:), sun3(:), sun4(:)};
sun_names = {'sun1', 'sun2', 'sun3', 'sun4'};

classicload = classicload(:);
fenshijijia = fenshijijia(:);
n = numel(wind_vectors{1});

% 目标：最大化 sum(t) = 最小化 -sum(t)，线性规划，使用 linprog
opts = optimoptions('linprog', 'Display', 'iter');
f  = -ones(n, 1);
lb = 0.1 * ones(n, 1);

problem4_results = struct([]);
result_idx = 0;

for wind_idx = 1:numel(wind_vectors)
    wind = wind_vectors{wind_idx};
    wind_name = wind_names{wind_idx};

    for sun_idx = 1:numel(sun_vectors)
        sun = sun_vectors{sun_idx};
        sun_name = sun_names{sun_idx};

        fprintf('\n==================== 情况 %02d：%s + %s ====================\n', ...
                (wind_idx - 1)*numel(sun_vectors) + sun_idx, wind_name, sun_name);

        % 约束 load <= generate：t*41.5 <= wind*40+sun*64 - classicload*6
        % 等价于收紧各时段上界
        ub = min(ones(n,1), (wind*40 + sun*64 - classicload*6) / (20.75*2));

        if any(ub < lb)
            warning('%s + %s：部分时段发电量不足以支撑最低制氢负荷，跳过。', wind_name, sun_name);
            continue;
        end

        t_opt = linprog(f, [], [], [], [], lb, ub, opts);
        t_opt = t_opt(:);
        sum_t = sum(t_opt);
        fprintf('最大 sum(t) = %.6g\n', sum_t);
        disp('最优 t 向量：'); disp(t_opt);

        metrics = calcMetrics(t_opt, classicload, wind, sun, fenshijijia);

        fprintf('新能源自发自用电量占总可用发电量比例 = %.6g\n', metrics.xnyzf);
        fprintf('总用电量中绿电占比 = %.6g\n', metrics.zydl);
        fprintf('新能源上网电量比例 = %.6g\n', metrics.xnysw);

        result_idx = result_idx + 1;
        problem4_results(result_idx).case_idx    = result_idx;
        problem4_results(result_idx).wind_idx    = wind_idx;
        problem4_results(result_idx).wind_name   = wind_name;
        problem4_results(result_idx).sun_idx     = sun_idx;
        problem4_results(result_idx).sun_name    = sun_name;
        problem4_results(result_idx).t_opt       = t_opt;
        problem4_results(result_idx).sum_t       = sum_t;
        problem4_results(result_idx).nh3load     = metrics.nh3load;
        problem4_results(result_idx).load        = metrics.load;
        problem4_results(result_idx).generate    = metrics.generate;
        problem4_results(result_idx).buy         = metrics.buy;
        problem4_results(result_idx).sell        = metrics.sell;
        problem4_results(result_idx).buycost     = metrics.buycost;
        problem4_results(result_idx).dayload     = metrics.dayload;
        problem4_results(result_idx).daygenerate = metrics.daygenerate;
        problem4_results(result_idx).daybuy      = metrics.daybuy;
        problem4_results(result_idx).daysell     = metrics.daysell;
        problem4_results(result_idx).xnyzf       = metrics.xnyzf;
        problem4_results(result_idx).zydl        = metrics.zydl;
        problem4_results(result_idx).xnysw       = metrics.xnysw;
    end
end

plotProblem4Distributions(problem4_results);

% -------------------------------------------------------------------------
function metrics = calcMetrics(t, classicload, wind, sun, fenshijijia)
    t = t(:); classicload = classicload(:);
    wind = wind(:); sun = sun(:); fenshijijia = fenshijijia(:);

    metrics.nh3load  = t * 20.75 * 2;
    metrics.load     = classicload * 6 + metrics.nh3load;
    metrics.generate = wind * 40 + sun * 64;

    diff = metrics.load - metrics.generate;
    metrics.buy  = diff; metrics.sell = diff;
    metrics.buy(metrics.buy < 0)   = 0;
    metrics.sell(metrics.sell > 0) = 0;

    metrics.buycost     = sum(metrics.buy .* fenshijijia);
    metrics.dayload     = sum(metrics.load);
    metrics.daygenerate = sum(metrics.generate);
    metrics.daybuy      = sum(metrics.buy);
    metrics.daysell     = abs(sum(metrics.sell));

    metrics.xnyzf = (metrics.dayload - metrics.daysell - metrics.daybuy) / metrics.daygenerate;
    metrics.zydl  = (metrics.daygenerate - metrics.daysell) / metrics.dayload;
    metrics.xnysw = metrics.daysell / metrics.daygenerate;
end

% -------------------------------------------------------------------------
function plotProblem4Distributions(problem4_results)
    if isempty(problem4_results)
        warning('没有可绘制的结果。');
        return;
    end

    sort_keys = [[problem4_results.wind_idx]', [problem4_results.sun_idx]'];
    [~, order] = sortrows(sort_keys, [1 2]);
    problem4_results = problem4_results(order);

    combo_labels = arrayfun(@(r) sprintf('%s+%s', r.wind_name, r.sun_name), ...
                            problem4_results, 'UniformOutput', false);

    metric_fields = {'sum_t', 'xnyzf', 'zydl', 'xnysw', 'daybuy', 'daysell'};
    metric_names  = {'sum(t)', 'problem4\_xnyzf', 'problem4\_zydl', ...
                     'problem4\_xnysw', 'problem4\_daybuy', 'problem4\_daysell'};

    figure('Name', 'problem4：各 wind/sun 组合指标分布', 'Color', 'w');

    for metric_idx = 1:numel(metric_fields)
        metric_field = metric_fields{metric_idx};
        metric_name  = metric_names{metric_idx};
        values = arrayfun(@(r) r.(metric_field), problem4_results);

        subplot(2, 3, metric_idx);
        bar(values, 'FaceColor', [0.2 0.45 0.75]);
        grid on; hold on;

        valid_values = values(~isnan(values));
        mean_value = mean(valid_values);
        min_value  = min(valid_values);
        max_value  = max(valid_values);
        [~, min_idx] = min(values);
        [~, max_idx] = max(values);

        plot(xlim, [mean_value mean_value], 'r--', 'LineWidth', 1.2);
        scatter(min_idx, min_value, 45, 'g', 'filled');
        scatter(max_idx, max_value, 45, 'r', 'filled');

        title({metric_name, sprintf('均值=%.4g，最小=%.4g，最大=%.4g', ...
               mean_value, min_value, max_value)}, 'Interpreter', 'none');
        xlabel('wind + sun 组合');
        ylabel(metric_name, 'Interpreter', 'none');
        set(gca, 'XTick', 1:numel(combo_labels), ...
                 'XTickLabel', combo_labels, 'XTickLabelRotation', 45);
        xlim([0.5, numel(combo_labels) + 0.5]);

        if min_value == max_value
            ylim([min_value - 0.5, max_value + 0.5]);
        end

        legend({'组合值', '均值', '最小值', '最大值'}, 'Location', 'best', 'FontSize', 8);
        hold off;
    end
end
