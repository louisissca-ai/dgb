wind_vectors = {wind1(:), wind2(:), wind3(:), wind4(:), wind5(:), wind6(:)};
wind_names = {'wind1', 'wind2', 'wind3', 'wind4', 'wind5', 'wind6'};
sun_vectors = {sun1(:), sun2(:), sun3(:), sun4(:)};
sun_names = {'sun1', 'sun2', 'sun3', 'sun4'};

classicload = classicload(:);
fenshijijia = fenshijijia(:);
n = numel(wind_vectors{1});


opts = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');


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




set(0, 'DefaultAxesFontName', 'SimHei');         
set(0, 'DefaultTextFontName', 'SimHei');
set(0, 'DefaultAxesFontSize', 10);

DAYS_PER_SCENARIO = 15;
N_SCENARIO        = 24;
N_DAY             = DAYS_PER_SCENARIO * N_SCENARIO;   % 360
beq_to_yield      = containers.Map([24 21 18 15 12], [72 63 54 45 36]);


TH_XNYZF = 0.60;
TH_ZYDL  = 0.30;
TH_XNYSW = 0.20;


n_beq = numel(beq_values);
result_grid = cell(N_SCENARIO, n_beq);

for k = 1:numel(problem3_results)
    r = problem3_results(k);
    scen_idx = (r.wind_idx - 1) * 4 + r.sun_idx;
    [~, beq_col] = ismember(r.beq, beq_values);
    result_grid{scen_idx, beq_col} = r;
end


J_mat      = zeros(N_SCENARIO, n_beq);   % 吨氨成本
xnyzf_mat  = zeros(N_SCENARIO, n_beq);
zydl_mat   = zeros(N_SCENARIO, n_beq);
xnysw_mat  = zeros(N_SCENARIO, n_beq);
buy_mat    = zeros(N_SCENARIO, n_beq);   % 日购电量
sell_mat   = zeros(N_SCENARIO, n_beq);   % 日上网电量
load_mat   = zeros(N_SCENARIO, n_beq);   % 日用电量
gen_mat    = zeros(N_SCENARIO, n_beq);   % 日发电量
scen_label = cell(N_SCENARIO, 1);

for s = 1:N_SCENARIO
    for b = 1:n_beq
        r = result_grid{s, b};
        J_mat(s, b)     = r.J_min;
        xnyzf_mat(s, b) = r.xnyzf;
        zydl_mat(s, b)  = r.zydl;
        xnysw_mat(s, b) = r.xnysw;
        buy_mat(s, b)   = r.daybuy;
        sell_mat(s, b)  = r.daysell;
        load_mat(s, b)  = r.dayload;
        gen_mat(s, b)   = r.daygenerate;
        if b == 1
            scen_label{s} = sprintf('%s+%s', r.wind_name, r.sun_name);
        end
    end
end

% 产量标签
yield_label = arrayfun(@(b) sprintf('%d吨/日', beq_to_yield(b)), beq_values, ...
                      'UniformOutput', false);

fig1 = figure('Name', '问题三(1) 全年吨氨成本分布曲线', ...
              'Color', 'w', 'Position', [100 100 1200 600]);

colors_5 = [0.85 0.33 0.10;    % 72 红
            0.93 0.69 0.13;    % 63 橙
            0.47 0.67 0.19;    % 54 绿
            0.30 0.75 0.93;    % 45 浅蓝
            0.00 0.45 0.74];   % 36 深蓝

day_axis = 1:N_DAY;
all_year_cost = zeros(N_DAY, n_beq);   % 360天 × 5产量

for b = 1:n_beq
    for s = 1:N_SCENARIO
        idx_start = (s - 1) * DAYS_PER_SCENARIO + 1;
        idx_end   = s * DAYS_PER_SCENARIO;
        all_year_cost(idx_start:idx_end, b) = J_mat(s, b);
    end
end


subplot(2, 1, 1);
hold on; grid on;
for b = 1:n_beq
    plot(day_axis, all_year_cost(:, b), '-', 'LineWidth', 1.5, ...
         'Color', colors_5(b, :), 'DisplayName', yield_label{b});
end

for s = 1:N_SCENARIO-1
    xline(s * DAYS_PER_SCENARIO + 0.5, ':', 'Color', [0.7 0.7 0.7], ...
          'HandleVisibility', 'off');
end
xlabel('全年天数（1~360）');
ylabel('吨氨成本 / (元·吨^{-1})');
title('图 1(a)  全年吨氨成本分布曲线（按场景顺序排列，每场景 15 天）');
legend('Location', 'eastoutside', 'NumColumns', 1);
xlim([0 N_DAY+1]);

subplot(2, 1, 2);
hold on; grid on;
for b = 1:n_beq
    sorted_cost = sort(all_year_cost(:, b), 'ascend');
    plot(day_axis, sorted_cost, '-', 'LineWidth', 1.8, ...
         'Color', colors_5(b, :), 'DisplayName', yield_label{b});
end
xlabel('累计天数（按吨氨成本升序排列）');
ylabel('吨氨成本 / (元·吨^{-1})');
title('图 1(b)  全年吨氨成本排序曲线（累计分布）');
legend('Location', 'eastoutside', 'NumColumns', 1);
xlim([0 N_DAY+1]);


saveas(fig1, 'problem3_annual_cost_curve.png');
fprintf('\n[图 1] 已保存：problem3_annual_cost_curve.png\n');


fig2 = figure('Name', '问题三(1) 吨氨成本箱线图', ...
              'Color', 'w', 'Position', [100 100 800 500]);


boxplot(J_mat, 'Labels', yield_label, 'Colors', colors_5, 'Widths', 0.6);
grid on;
xlabel('日产量方案');
ylabel('吨氨成本 / (元·吨^{-1})');
title('图 2  24 风光场景下不同产量方案的吨氨成本分布');

hold on;
mean_J = mean(J_mat, 1);
plot(1:n_beq, mean_J, 'd', 'MarkerSize', 8, ...
     'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
for b = 1:n_beq
    text(b, mean_J(b), sprintf('  μ=%.0f', mean_J(b)), ...
         'FontSize', 9, 'Color', 'k');
end
saveas(fig2, 'problem3_cost_boxplot.png');
fprintf('[图 2] 已保存：problem3_cost_boxplot.png\n');


fig3 = figure('Name', '问题三(1) 绿电指标合规热图', ...
              'Color', 'w', 'Position', [100 100 1400 500]);

pass_xnyzf = xnyzf_mat > TH_XNYZF;        % 自发自用 > 60%
pass_zydl  = zydl_mat  > TH_ZYDL;          % 绿电占比 > 30%
pass_xnysw = xnysw_mat < TH_XNYSW;        % 上网比例 < 20%



subplot(1, 3, 1);
imagesc(double(pass_xnyzf));
colormap(gca, [0.85 0.33 0.10; 0.47 0.67 0.19]);    
caxis([0 1]);
title({'自发自用电量占比 (>60%)', sprintf('达标率 %.1f%%', 100*mean(pass_xnyzf(:)))});
xlabel('日产量方案'); ylabel('风光场景');
set(gca, 'XTick', 1:n_beq, 'XTickLabel', yield_label, ...
         'YTick', 1:N_SCENARIO, 'YTickLabel', scen_label, 'FontSize', 8);
xtickangle(45);

subplot(1, 3, 2);
imagesc(double(pass_zydl));
colormap(gca, [0.85 0.33 0.10; 0.47 0.67 0.19]);
caxis([0 1]);
title({'总用电量绿电占比 (>30%)', sprintf('达标率 %.1f%%', 100*mean(pass_zydl(:)))});
xlabel('日产量方案'); ylabel('风光场景');
set(gca, 'XTick', 1:n_beq, 'XTickLabel', yield_label, ...
         'YTick', 1:N_SCENARIO, 'YTickLabel', scen_label, 'FontSize', 8);
xtickangle(45);

subplot(1, 3, 3);
imagesc(double(pass_xnysw));
colormap(gca, [0.85 0.33 0.10; 0.47 0.67 0.19]);
caxis([0 1]);
title({'新能源上网比例 (<20%)', sprintf('达标率 %.1f%%', 100*mean(pass_xnysw(:)))});
xlabel('日产量方案'); ylabel('风光场景');
set(gca, 'XTick', 1:n_beq, 'XTickLabel', yield_label, ...
         'YTick', 1:N_SCENARIO, 'YTickLabel', scen_label, 'FontSize', 8);
xtickangle(45);

sgtitle('图 3  绿电直连三指标合规情况（绿色=达标，红色=不达标）');
saveas(fig3, 'problem3_compliance_heatmap.png');
fprintf('[图 3] 已保存：problem3_compliance_heatmap.png\n');


fig4 = figure('Name', '问题三(1) 三类合规统计', ...
              'Color', 'w', 'Position', [100 100 1000 500]);

cat_all_pass   = (pass_count == 3);
cat_partial    = (pass_count >= 1) & (pass_count <= 2);
cat_none_pass  = (pass_count == 0);


days_all    = sum(cat_all_pass,  1) * DAYS_PER_SCENARIO;
days_part   = sum(cat_partial,   1) * DAYS_PER_SCENARIO;
days_none   = sum(cat_none_pass, 1) * DAYS_PER_SCENARIO;

subplot(1, 2, 1);
bar_data = [days_all; days_part; days_none]';   
b_h = bar(bar_data, 'stacked');
b_h(1).FaceColor = [0.47 0.67 0.19];  
b_h(2).FaceColor = [0.93 0.69 0.13];   
b_h(3).FaceColor = [0.85 0.33 0.10];   
grid on;
xlabel('日产量方案'); ylabel('全年天数');
set(gca, 'XTickLabel', yield_label);
legend({'全满足', '部分满足', '全不满足'}, 'Location', 'best');
title('图 4(a)  各产量方案全年绿电指标合规天数');
ylim([0 N_DAY*1.05]);
% 在柱顶标注
for i = 1:n_beq
    text(i, N_DAY*1.02, sprintf('%d/%d/%d', days_all(i), days_part(i), days_none(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 9);
end

subplot(1, 2, 2);
total_all  = sum(days_all);
total_part = sum(days_part);
total_none = sum(days_none);
pie_data = [total_all, total_part, total_none];
pie_labels = {sprintf('全满足 %d天 (%.1f%%)', total_all,  100*total_all /sum(pie_data)), ...
              sprintf('部分满足 %d天 (%.1f%%)', total_part, 100*total_part/sum(pie_data)), ...
              sprintf('全不满足 %d天 (%.1f%%)', total_none, 100*total_none/sum(pie_data))};
pie_h = pie(pie_data, pie_labels);
colormap_local = [0.47 0.67 0.19; 0.93 0.69 0.13; 0.85 0.33 0.10];
for i = 1:3
    pie_h(2*i-1).FaceColor = colormap_local(i, :);
end
title('图 4(b)  全部产量方案合并的全年合规分布');

saveas(fig4, 'problem3_compliance_summary.png');
fprintf('[图 4] 已保存：problem3_compliance_summary.png\n');


fig5 = figure('Name', '问题三(2) 购售电量分布', ...
              'Color', 'w', 'Position', [100 100 1200 500]);

subplot(1, 2, 1);
boxplot(buy_mat, 'Labels', yield_label, 'Widths', 0.6);
grid on;
xlabel('日产量方案'); ylabel('日购电量 / (MW·h)');
title('图 5(a)  24 场景下日购电量分布');
hold on;
plot(1:n_beq, mean(buy_mat, 1), 'rd', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

subplot(1, 2, 2);
boxplot(sell_mat, 'Labels', yield_label, 'Widths', 0.6);
grid on;
xlabel('日产量方案'); ylabel('日上网电量 / (MW·h)');
title('图 5(b)  24 场景下日上网电量分布');
hold on;
plot(1:n_beq, mean(sell_mat, 1), 'rd', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

saveas(fig5, 'problem3_buy_sell_distribution.png');
fprintf('[图 5] 已保存：problem3_buy_sell_distribution.png\n');


fig6 = figure('Name', '问题三(2) 绿电指标分布', ...
              'Color', 'w', 'Position', [100 100 1500 450]);

subplot(1, 3, 1);
boxplot(xnyzf_mat, 'Labels', yield_label, 'Widths', 0.6);
yline(TH_XNYZF, 'r--', sprintf('阈值 %.0f%%', TH_XNYZF*100), 'LineWidth', 1.5);
grid on; ylabel('自发自用占比');
title('图 6(a) 新能源自发自用电量占比');
ylim([0 1.05]);

subplot(1, 3, 2);
boxplot(zydl_mat, 'Labels', yield_label, 'Widths', 0.6);
yline(TH_ZYDL, 'r--', sprintf('阈值 %.0f%%', TH_ZYDL*100), 'LineWidth', 1.5);
grid on; ylabel('绿电占比');
title('图 6(b) 总用电量绿电占比');

subplot(1, 3, 3);
boxplot(xnysw_mat, 'Labels', yield_label, 'Widths', 0.6);
yline(TH_XNYSW, 'r--', sprintf('阈值 %.0f%%', TH_XNYSW*100), 'LineWidth', 1.5);
grid on; ylabel('上网比例');
title('图 6(c) 新能源上网电量比例');

saveas(fig6, 'problem3_three_indicators.png');
fprintf('[图 6] 已保存：problem3_three_indicators.png\n');

fprintf('\n\n');
fprintf('=========================================================\n');
fprintf('              问题三 关键数值结论汇总                    \n');
fprintf('=========================================================\n');

fprintf('\n【1】各产量方案吨氨成本统计（元/吨）：\n');
fprintf('  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s\n', ...
        '产量', '最小', '最大', '均值', '中位数', '标准差');
for b = 1:n_beq
    fprintf('  %-10s  %-10.2f  %-10.2f  %-10.2f  %-10.2f  %-10.2f\n', ...
            yield_label{b}, min(J_mat(:,b)), max(J_mat(:,b)), ...
            mean(J_mat(:,b)), median(J_mat(:,b)), std(J_mat(:,b)));
end


fprintf('\n【2】全年总吨氨成本（按产量方案分别核算）：\n');
fprintf('  %-10s  %-15s  %-15s  %-15s\n', ...
        '产量', '全年总成本(元)', '全年总产氨(吨)', '全年吨氨成本(元/吨)');
total_cost_per_yield = zeros(1, n_beq);
total_ammonia_per_yield = zeros(1, n_beq);
annual_cost_per_yield = zeros(1, n_beq);
for b = 1:n_beq
    yield_b = beq_to_yield(beq_values(b));    
    total_cost_per_yield(b)    = sum(J_mat(:, b) * yield_b * DAYS_PER_SCENARIO);
    total_ammonia_per_yield(b) = sum(yield_b * DAYS_PER_SCENARIO * ones(N_SCENARIO, 1));
    annual_cost_per_yield(b)   = total_cost_per_yield(b) / total_ammonia_per_yield(b);
    fprintf('  %-10s  %-15.2e  %-15.0f  %-15.2f\n', ...
            yield_label{b}, total_cost_per_yield(b), ...
            total_ammonia_per_yield(b), annual_cost_per_yield(b));
end

[best_annual_cost, best_b] = min(annual_cost_per_yield);
fprintf('\n  ==> 全年综合最优产量方案：%s，全年吨氨成本 %.2f 元/吨\n', ...
        yield_label{best_b}, best_annual_cost);


fprintf('\n【3】绿电指标全年合规分布：\n');
fprintf('  %-10s  %-12s  %-12s  %-12s\n', ...
        '产量', '全满足(天)', '部分满足(天)', '全不满足(天)');
for b = 1:n_beq
    fprintf('  %-10s  %-12d  %-12d  %-12d\n', ...
            yield_label{b}, days_all(b), days_part(b), days_none(b));
end
fprintf('  %-10s  %-12d  %-12d  %-12d\n', ...
        '合计', sum(days_all), sum(days_part), sum(days_none));


fprintf('\n【4】单项指标达标情况（按产量）：\n');
fprintf('  %-10s  %-15s  %-15s  %-15s\n', ...
        '产量', '自发自用达标率', '绿电占比达标率', '上网比例达标率');
for b = 1:n_beq
    fprintf('  %-10s  %-15.1f%%  %-15.1f%%  %-15.1f%%\n', ...
            yield_label{b}, ...
            100*mean(pass_xnyzf(:, b)), ...
            100*mean(pass_zydl(:, b)), ...
            100*mean(pass_xnysw(:, b)));
end


fprintf('\n【5】日购售电量统计（MW·h）：\n');
fprintf('  %-10s  %-12s  %-12s  %-12s  %-12s\n', ...
        '产量', '购电均值', '购电最大', '上网均值', '上网最大');
for b = 1:n_beq
    fprintf('  %-10s  %-12.2f  %-12.2f  %-12.2f  %-12.2f\n', ...
            yield_label{b}, ...
            mean(buy_mat(:, b)), max(buy_mat(:, b)), ...
            mean(sell_mat(:, b)), max(sell_mat(:, b)));
end


fprintf('\n【6】各产量方案下的最优/最差场景：\n');
fprintf('  %-10s  %-25s  %-25s\n', '产量', '最低成本场景', '最高成本场景');
for b = 1:n_beq
    [min_j, idx_min] = min(J_mat(:, b));
    [max_j, idx_max] = max(J_mat(:, b));
    fprintf('  %-10s  %-15s(%.0f)  %-15s(%.0f)\n', ...
            yield_label{b}, scen_label{idx_min}, min_j, ...
            scen_label{idx_max}, max_j);
end


fprintf('\n【7】综合建议：\n');
[~, best_mean_b] = min(mean(J_mat, 1));
fprintf('   - 平均吨氨成本最低的产量方案：%s（均值 %.2f 元/吨）\n', ...
        yield_label{best_mean_b}, mean(J_mat(:, best_mean_b)));
fprintf('   - 全年综合最优产量方案（按全年总成本核算）：%s（全年 %.2f 元/吨）\n', ...
        yield_label{best_b}, best_annual_cost);


if all(days_none > 0) || all(days_part > 0)
    fprintf('   - ⚠️ 所有产量方案均存在绿电指标不完全合规的天数，\n');
    fprintf('     说明在仅依赖风光直连+电网交互的运行模式下，\n');
    fprintf('     难以在全部场景实现"自发自用>60%% + 上网<20%%"双重约束。\n');
end

fprintf('\n=========================================================\n');
fprintf('  绘图与统计完成。所有图片已保存至当前工作目录。\n');
fprintf('=========================================================\n');


try
   
    out_J = array2table(J_mat, 'VariableNames', ...
        {'Y72', 'Y63', 'Y54', 'Y45', 'Y36'}, ...
        'RowNames', scen_label);
    writetable(out_J, 'problem3_summary.xlsx', 'Sheet', '吨氨成本', ...
               'WriteRowNames', true);

    out_xnyzf = array2table(xnyzf_mat, 'VariableNames', ...
        {'Y72', 'Y63', 'Y54', 'Y45', 'Y36'}, 'RowNames', scen_label);
    writetable(out_xnyzf, 'problem3_summary.xlsx', 'Sheet', '自发自用', ...
               'WriteRowNames', true);

    out_zydl = array2table(zydl_mat, 'VariableNames', ...
        {'Y72', 'Y63', 'Y54', 'Y45', 'Y36'}, 'RowNames', scen_label);
    writetable(out_zydl, 'problem3_summary.xlsx', 'Sheet', '绿电占比', ...
               'WriteRowNames', true);

    out_xnysw = array2table(xnysw_mat, 'VariableNames', ...
        {'Y72', 'Y63', 'Y54', 'Y45', 'Y36'}, 'RowNames', scen_label);
    writetable(out_xnysw, 'problem3_summary.xlsx', 'Sheet', '上网比例', ...
               'WriteRowNames', true);

    fprintf('\n数据已导出：problem3_summary.xlsx\n');
catch ME
    warning('Excel 导出失败：%s（图片均已成功保存）', ME.message);
end



set(0, 'DefaultAxesFontName', 'SimHei');
set(0, 'DefaultTextFontName', 'SimHei');
set(0, 'DefaultAxesFontSize', 10);

DAYS_PER_SCENARIO = 15;
N_SCENARIO        = 24;
N_DAY             = DAYS_PER_SCENARIO * N_SCENARIO;
beq_values_local  = [24, 21, 18, 15, 12];
beq_to_yield      = containers.Map([24 21 18 15 12], [72 63 54 45 36]);
n_beq             = numel(beq_values_local);

% 绿电三指标阈值
TH_XNYZF = 0.60;
TH_ZYDL  = 0.30;
TH_XNYSW = 0.20;

yield_label = arrayfun(@(b) sprintf('%d吨/日', beq_to_yield(b)), beq_values_local, ...
                      'UniformOutput', false);

[J2,  xnyzf2, zydl2, xnysw2, buy2, sell2, scen_label] = extractMat(problem2_results, beq_values_local);
[J3,  xnyzf3, zydl3, xnysw3, buy3, sell3, ~]          = extractMat(problem3_results, beq_values_local);

dJ      = J3      - J2;      
dXnyzf  = xnyzf3  - xnyzf2;
dZydl   = zydl3   - zydl2;
dXnysw  = xnysw3  - xnysw2;
dBuy    = buy3    - buy2;
dSell   = sell3   - sell2;


fig7 = figure('Name', '问题三(3) 吨氨成本对比', ...
              'Color', 'w', 'Position', [100 100 1200 500]);


subplot(1, 2, 1);
mean_J2 = mean(J2, 1);
mean_J3 = mean(J3, 1);
std_J2  = std(J2, 0, 1);
std_J3  = std(J3, 0, 1);

bar_data = [mean_J2; mean_J3]';   % 5×2
b_h = bar(bar_data, 'grouped');
b_h(1).FaceColor = [0.85 0.33 0.10];   % 问题二 红
b_h(2).FaceColor = [0.00 0.45 0.74];   % 问题三 蓝
hold on; grid on;


ngroups = n_beq;
nbars   = 2;
groupwidth = min(0.8, nbars/(nbars+1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    if i == 1
        errorbar(x, mean_J2, std_J2, 'k.', 'LineWidth', 1, 'CapSize', 8);
    else
        errorbar(x, mean_J3, std_J3, 'k.', 'LineWidth', 1, 'CapSize', 8);
    end
end

set(gca, 'XTickLabel', yield_label);
xlabel('日产量方案'); ylabel('吨氨成本均值 / (元·吨^{-1})');
legend({'问题二（离散开停）', '问题三（连续可调）'}, 'Location', 'best');
title('图 7(a)  吨氨成本均值对比（误差棒为24场景标准差）');


for i = 1:n_beq
    drop_pct = (mean_J2(i) - mean_J3(i)) / mean_J2(i) * 100;
    text(i, max(mean_J2(i), mean_J3(i)) + max(std_J2(i), std_J3(i)) + 50, ...
         sprintf('↓%.1f%%', drop_pct), ...
         'HorizontalAlignment', 'center', 'FontSize', 10, ...
         'Color', [0.47 0.67 0.19], 'FontWeight', 'bold');
end


subplot(1, 2, 2);
hold on; grid on;
colors_5 = [0.85 0.33 0.10; 0.93 0.69 0.13; 0.47 0.67 0.19; ...
            0.30 0.75 0.93; 0.00 0.45 0.74];
for b = 1:n_beq
    scatter(J2(:, b), J3(:, b), 50, colors_5(b, :), 'filled', ...
            'DisplayName', yield_label{b});
end


all_J = [J2(:); J3(:)];
lim_min = min(all_J) * 0.95;
lim_max = max(all_J) * 1.05;
plot([lim_min lim_max], [lim_min lim_max], 'k--', 'LineWidth', 1, ...
     'DisplayName', 'y = x（无差异线）');

xlabel('问题二 吨氨成本 / (元·吨^{-1})');
ylabel('问题三 吨氨成本 / (元·吨^{-1})');
title('图 7(b)  24场景逐点对比（位于y=x下方表示连续可调更优）');
legend('Location', 'best');
xlim([lim_min lim_max]); ylim([lim_min lim_max]);
axis square;

saveas(fig7, 'problem2vs3_cost_compare.png');
fprintf('[图 7] 已保存：problem2vs3_cost_compare.png\n');


fig8 = figure('Name', '问题三(3) 绿电指标对比', ...
              'Color', 'w', 'Position', [100 100 1500 500]);

plot_three_indicator_compare(1, xnyzf2, xnyzf3, yield_label, TH_XNYZF, ...
                              '图 8(a) 新能源自发自用占比', '自发自用占比');
plot_three_indicator_compare(2, zydl2, zydl3, yield_label, TH_ZYDL, ...
                              '图 8(b) 总用电绿电占比', '绿电占比');
plot_three_indicator_compare(3, xnysw2, xnysw3, yield_label, TH_XNYSW, ...
                              '图 8(c) 新能源上网比例', '上网比例');

saveas(fig8, 'problem2vs3_indicators_compare.png');
fprintf('[图 8] 已保存：problem2vs3_indicators_compare.png\n');

fig9 = figure('Name', '问题三(3) 合规天数对比', ...
              'Color', 'w', 'Position', [100 100 1200 500]);


[days_all_p2, days_part_p2, days_none_p2] = countCompliance( ...
    xnyzf2, zydl2, xnysw2, TH_XNYZF, TH_ZYDL, TH_XNYSW, DAYS_PER_SCENARIO);
[days_all_p3, days_part_p3, days_none_p3] = countCompliance( ...
    xnyzf3, zydl3, xnysw3, TH_XNYZF, TH_ZYDL, TH_XNYSW, DAYS_PER_SCENARIO);


subplot(1, 2, 1);
bar_data = [days_all_p2; days_part_p2; days_none_p2]';
b_h = bar(bar_data, 'stacked');
b_h(1).FaceColor = [0.47 0.67 0.19];
b_h(2).FaceColor = [0.93 0.69 0.13];
b_h(3).FaceColor = [0.85 0.33 0.10];
grid on;
xlabel('日产量方案'); ylabel('全年天数');
set(gca, 'XTickLabel', yield_label);
legend({'全满足', '部分满足', '全不满足'}, 'Location', 'best');
title('图 9(a)  问题二（离散开停）合规分布');
ylim([0 N_DAY*1.05]);


subplot(1, 2, 2);
bar_data = [days_all_p3; days_part_p3; days_none_p3]';
b_h = bar(bar_data, 'stacked');
b_h(1).FaceColor = [0.47 0.67 0.19];
b_h(2).FaceColor = [0.93 0.69 0.13];
b_h(3).FaceColor = [0.85 0.33 0.10];
grid on;
xlabel('日产量方案'); ylabel('全年天数');
set(gca, 'XTickLabel', yield_label);
legend({'全满足', '部分满足', '全不满足'}, 'Location', 'best');
title('图 9(b)  问题三（连续可调）合规分布');
ylim([0 N_DAY*1.05]);

saveas(fig9, 'problem2vs3_compliance_compare.png');
fprintf('[图 9] 已保存：problem2vs3_compliance_compare.png\n');


fig10 = figure('Name', '问题三(3) 成本下降热图', ...
              'Color', 'w', 'Position', [100 100 800 700]);


dJ_pct = (J2 - J3) ./ J2 * 100;  

imagesc(dJ_pct);
colormap(redblueColormap());

caxis_max = max(abs(dJ_pct(:)));
caxis([-caxis_max, caxis_max]);
cb = colorbar;
cb.Label.String = '吨氨成本下降百分比 (%)';
cb.Label.FontSize = 11;

for s = 1:N_SCENARIO
    for b = 1:n_beq
        text(b, s, sprintf('%.1f', dJ_pct(s, b)), ...
             'HorizontalAlignment', 'center', 'FontSize', 7, ...
             'Color', 'k');
    end
end

xlabel('日产量方案'); ylabel('风光场景');
set(gca, 'XTick', 1:n_beq, 'XTickLabel', yield_label, ...
         'YTick', 1:N_SCENARIO, 'YTickLabel', scen_label, 'FontSize', 8);
xtickangle(45);
title({'图 10  连续可调相对离散开停的吨氨成本下降幅度 (%)', ...
       '蓝色：连续可调更优；红色：离散开停更优'});

saveas(fig10, 'problem2vs3_cost_drop_heatmap.png');
fprintf('[图 10] 已保存：problem2vs3_cost_drop_heatmap.png\n');

fig11 = figure('Name', '问题三(3) 全年成本曲线对比', ...
              'Color', 'w', 'Position', [100 100 1200 700]);

day_axis = 1:N_DAY;

for b = 1:n_beq
    subplot(n_beq, 1, b);
    hold on; grid on;
    
    cost_p2 = zeros(N_DAY, 1);
    cost_p3 = zeros(N_DAY, 1);
    for s = 1:N_SCENARIO
        idx_start = (s - 1) * DAYS_PER_SCENARIO + 1;
        idx_end   = s * DAYS_PER_SCENARIO;
        cost_p2(idx_start:idx_end) = J2(s, b);
        cost_p3(idx_start:idx_end) = J3(s, b);
    end
    
    plot(day_axis, cost_p2, '-', 'LineWidth', 1.2, ...
         'Color', [0.85 0.33 0.10], 'DisplayName', '离散开停（问题二）');
    plot(day_axis, cost_p3, '-', 'LineWidth', 1.2, ...
         'Color', [0.00 0.45 0.74], 'DisplayName', '连续可调（问题三）');
    
    for s = 1:N_SCENARIO-1
        xline(s * DAYS_PER_SCENARIO + 0.5, ':', 'Color', [0.85 0.85 0.85], ...
              'HandleVisibility', 'off');
    end
    
    ylabel(sprintf('%s\n吨氨成本(元/吨)', yield_label{b}));
    if b == 1
        title('图 11  两种模式全年吨氨成本曲线叠加对比（5种产量分行）');
        legend('Location', 'best');
    end
    if b == n_beq
        xlabel('全年天数（1~360）');
    else
        set(gca, 'XTickLabel', []);
    end
    xlim([0 N_DAY+1]);
end

saveas(fig11, 'problem2vs3_annual_curve_compare.png');
fprintf('[图 11] 已保存：problem2vs3_annual_curve_compare.png\n');


fprintf('\n\n');
fprintf('=========================================================\n');
fprintf('       问题三(3) 连续可调 vs 离散开停 对比结论          \n');
fprintf('=========================================================\n');

fprintf('\n【1】吨氨成本对比（按产量分组，单位：元/吨）：\n');
fprintf('  %-10s  %-12s  %-12s  %-12s  %-10s\n', ...
        '产量', '问题二均值', '问题三均值', '下降值', '下降幅度');
for b = 1:n_beq
    drop_abs = mean_J2(b) - mean_J3(b);
    drop_pct = drop_abs / mean_J2(b) * 100;
    fprintf('  %-10s  %-12.2f  %-12.2f  %-12.2f  %-10.2f%%\n', ...
            yield_label{b}, mean_J2(b), mean_J3(b), drop_abs, drop_pct);
end

fprintf('\n【2】全年总吨氨成本对比：\n');
fprintf('  %-10s  %-15s  %-15s  %-12s\n', ...
        '产量', '问题二(元/吨)', '问题三(元/吨)', '下降幅度');
for b = 1:n_beq
    yield_b = beq_to_yield(beq_values_local(b));
    tot_p2 = sum(J2(:, b) * yield_b * DAYS_PER_SCENARIO) / ...
             (yield_b * DAYS_PER_SCENARIO * N_SCENARIO);
    tot_p3 = sum(J3(:, b) * yield_b * DAYS_PER_SCENARIO) / ...
             (yield_b * DAYS_PER_SCENARIO * N_SCENARIO);
    fprintf('  %-10s  %-15.2f  %-15.2f  %-12.2f%%\n', ...
            yield_label{b}, tot_p2, tot_p3, (tot_p2 - tot_p3) / tot_p2 * 100);
end

fprintf('\n【3】绿电三指标均值变化（连续 - 离散）：\n');
fprintf('  %-10s  %-15s  %-15s  %-15s\n', ...
        '产量', 'Δ自发自用', 'Δ绿电占比', 'Δ上网比例');
for b = 1:n_beq
    fprintf('  %-10s  %-15.4f  %-15.4f  %-15.4f\n', ...
            yield_label{b}, ...
            mean(dXnyzf(:, b)), mean(dZydl(:, b)), mean(dXnysw(:, b)));
end

fprintf('\n【4】全年合规天数对比（全满足/部分满足/全不满足）：\n');
fprintf('  %-10s  %-25s  %-25s\n', '产量', '问题二（离散）', '问题三（连续）');
for b = 1:n_beq
    fprintf('  %-10s  %3d / %3d / %3d %10s %3d / %3d / %3d\n', ...
            yield_label{b}, ...
            days_all_p2(b), days_part_p2(b), days_none_p2(b), '', ...
            days_all_p3(b), days_part_p3(b), days_none_p3(b));
end
fprintf('  %-10s  %3d / %3d / %3d %10s %3d / %3d / %3d\n', ...
        '合计', ...
        sum(days_all_p2), sum(days_part_p2), sum(days_none_p2), '', ...
        sum(days_all_p3), sum(days_part_p3), sum(days_none_p3));

fprintf('\n【5】日购售电量变化（均值，连续 - 离散）：\n');
fprintf('  %-10s  %-15s  %-15s\n', '产量', 'Δ购电(MW·h)', 'Δ上网(MW·h)');
for b = 1:n_beq
    fprintf('  %-10s  %-15.3f  %-15.3f\n', ...
            yield_label{b}, mean(dBuy(:, b)), mean(dSell(:, b)));
end

fprintf('\n【6】综合结论：\n');
total_drop_pct = (mean(J2(:)) - mean(J3(:))) / mean(J2(:)) * 100;
fprintf('   - 连续可调相对离散开停，整体吨氨成本平均下降 %.2f%%。\n', total_drop_pct);

[~, best_drop_b] = max(mean_J2 - mean_J3);
fprintf('   - 受益最大的产量方案：%s（绝对下降 %.2f 元/吨）。\n', ...
        yield_label{best_drop_b}, mean_J2(best_drop_b) - mean_J3(best_drop_b));

compliance_imp_all  = sum(days_all_p3)  - sum(days_all_p2);
compliance_imp_none = sum(days_none_p2) - sum(days_none_p3);
fprintf('   - 全年"全满足"天数增加 %d 天，"全不满足"天数减少 %d 天。\n', ...
        compliance_imp_all, compliance_imp_none);

if mean(dSell(:)) < 0
    fprintf('   - 上网电量平均减少 %.2f MW·h/日，说明连续调节提升了新能源就地消纳能力。\n', ...
            -mean(dSell(:)));
end
if mean(dBuy(:)) < 0
    fprintf('   - 购电量平均减少 %.2f MW·h/日，说明连续调节优化了电网交互。\n', ...
            -mean(dBuy(:)));
end

fprintf('\n=========================================================\n');
fprintf('  问题二 vs 问题三 对比完成。\n');
fprintf('=========================================================\n');

function [J_mat, xnyzf_mat, zydl_mat, xnysw_mat, buy_mat, sell_mat, scen_label] = ...
         extractMat(results, beq_values_local)
    N_SCENARIO = 24;
    n_beq = numel(beq_values_local);
    J_mat     = zeros(N_SCENARIO, n_beq);
    xnyzf_mat = zeros(N_SCENARIO, n_beq);
    zydl_mat  = zeros(N_SCENARIO, n_beq);
    xnysw_mat = zeros(N_SCENARIO, n_beq);
    buy_mat   = zeros(N_SCENARIO, n_beq);
    sell_mat  = zeros(N_SCENARIO, n_beq);
    scen_label = cell(N_SCENARIO, 1);

    for k = 1:numel(results)
        r = results(k);
        scen_idx = (r.wind_idx - 1) * 4 + r.sun_idx;
        [~, beq_col] = ismember(r.beq, beq_values_local);
        J_mat(scen_idx, beq_col)     = r.J_min;
        xnyzf_mat(scen_idx, beq_col) = r.xnyzf;
        zydl_mat(scen_idx, beq_col)  = r.zydl;
        xnysw_mat(scen_idx, beq_col) = r.xnysw;
        buy_mat(scen_idx, beq_col)   = r.daybuy;
        sell_mat(scen_idx, beq_col)  = r.daysell;
        if beq_col == 1
            scen_label{scen_idx} = sprintf('%s+%s', r.wind_name, r.sun_name);
        end
    end
end

function plot_three_indicator_compare(sub_idx, mat2, mat3, yield_label, threshold, ...
                                      ttl, ylbl)
    subplot(1, 3, sub_idx);
    n_beq = size(mat2, 2);
    plot_data = [];
    plot_group = [];
    color_grp = [];
    for b = 1:n_beq
        plot_data  = [plot_data;  mat2(:, b); mat3(:, b)];
        plot_group = [plot_group; (2*b-1)*ones(size(mat2,1),1); 2*b*ones(size(mat3,1),1)];
        color_grp  = [color_grp;  ones(size(mat2,1),1);          2*ones(size(mat3,1),1)];
    end
    
    positions = 1:(2*n_beq);
    boxplot(plot_data, plot_group, 'Positions', positions, 'Widths', 0.6, ...
            'Colors', repmat([0.85 0.33 0.10; 0.00 0.45 0.74], n_beq, 1));
    
    yline(threshold, 'r--', sprintf('阈值 %.0f%%', threshold*100), ...
          'LineWidth', 1.5);
    grid on;
    
    tick_pos = 1.5:2:(2*n_beq);
    set(gca, 'XTick', tick_pos, 'XTickLabel', yield_label);
    xlabel('日产量方案（每对：左=离散，右=连续）');
    ylabel(ylbl);
    title(ttl);
end

function [days_all, days_part, days_none] = countCompliance( ...
         xnyzf_mat, zydl_mat, xnysw_mat, TH_XNYZF, TH_ZYDL, TH_XNYSW, ...
         DAYS_PER_SCENARIO)
    pass_xnyzf = xnyzf_mat > TH_XNYZF;
    pass_zydl  = zydl_mat  > TH_ZYDL;
    pass_xnysw = xnysw_mat < TH_XNYSW;
    pass_count = double(pass_xnyzf) + double(pass_zydl) + double(pass_xnysw);
    
    cat_all_pass   = (pass_count == 3);
    cat_partial    = (pass_count >= 1) & (pass_count <= 2);
    cat_none_pass  = (pass_count == 0);
    
    days_all  = sum(cat_all_pass,  1) * DAYS_PER_SCENARIO;
    days_part = sum(cat_partial,   1) * DAYS_PER_SCENARIO;
    days_none = sum(cat_none_pass, 1) * DAYS_PER_SCENARIO;
end

function cmap = redblueColormap()
    n = 256;
    half = n / 2;
    r = [linspace(0.7, 1, half), linspace(1, 0.0, half)];
    g = [linspace(0.0, 1, half), linspace(1, 0.45, half)];
    b = [linspace(0.0, 1, half), linspace(1, 0.74, half)];
    cmap = [r', g', b'];
end
