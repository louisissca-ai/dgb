wind_vectors = {wind1(:), wind2(:), wind3(:), wind4(:), wind5(:), wind6(:)};
wind_names = {'wind1', 'wind2', 'wind3', 'wind4', 'wind5', 'wind6'};
sun_vectors = {sun1(:), sun2(:), sun3(:), sun4(:)};
sun_names = {'sun1', 'sun2', 'sun3', 'sun4'};

classicload = classicload(:);
n = numel(wind_vectors{1});

% 目标：最大化 sum(t) = 最小化 -sum(t)，线性规划，使用 linprog
opts = optimoptions('linprog', 'Display', 'iter');
f  = -ones(n, 1);
lb = 0.1 * ones(n, 1);

problem4_results = struct([]);
result_idx = 0;

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
            if ub_val(i) < 0.1
                slack(i) = ub_val(i) - 0.1;   % 缺电：负值
            elseif raw_ub(i) > 1
                slack(i) = raw_ub(i) - 1;      % 弃电：正值
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


        J = objFun(t_opt, wind, sun);
        fprintf('目标函数 J = %.6g\n', J);

        result_idx = result_idx + 1;
        problem4_results(result_idx).case_idx    = result_idx;
        problem4_results(result_idx).wind_idx    = wind_idx;
        problem4_results(result_idx).wind_name   = wind_name;
        problem4_results(result_idx).sun_idx     = sun_idx;
        problem4_results(result_idx).sun_name    = sun_name;
        problem4_results(result_idx).t_opt       = t_opt;
        problem4_results(result_idx).sum_t       = sum_t;
        % ---------------- 新增：年产量与利用率指标 ----------------
        % 制氢氨满载等效功率 (MW)：碱性10 + PEM10 + 合成氨0.75 = 20.75，产能扩容2倍 → 41.5
        P_h2_total      = 20.75 * 2;
        annual_nh3      = sum_t * 3 * 365;                    % 全年制氨总量 (吨)
                                                              %   3 = 1.5吨/h × 1h步长 × 2倍扩容
        cap_utilization = sum_t / n;                          % 氢氨年平均产能利用率 = mean(t)
        gen_sum         = sum(wind*40 + sun*64);              % 日发电总量 (MW·h)
        power_util      = sum_t * P_h2_total / gen_sum;       % 发电利用率：制氨耗电 / 风光发电

        problem4_results(result_idx).annual_nh3      = annual_nh3;
        problem4_results(result_idx).cap_utilization = cap_utilization;
        problem4_results(result_idx).power_util      = power_util;
        problem4_results(result_idx).gen_sum         = gen_sum;

        fprintf('  全年制氨总量    = %.2f 吨\n', annual_nh3);
        fprintf('  产能利用率      = %.4f\n',    cap_utilization);
        fprintf('  发电利用率      = %.4f\n',    power_util);
    end
end



sort_keys = [[problem4_results.wind_idx]', [problem4_results.sun_idx]'];
[~, p1_order] = sortrows(sort_keys, [1 2]);
p1_sorted = problem4_results(p1_order);
combo_labels_p1 = arrayfun(@(r) sprintf('%s+%s', r.wind_name, r.sun_name), ...
                           p1_sorted, 'UniformOutput', false);

annual_vals = [p1_sorted.annual_nh3];
cap_vals    = [p1_sorted.cap_utilization];
power_vals  = [p1_sorted.power_util];

figure('Name', 'problem4-1：年产量与利用率分布', 'Color', 'w', ...
       'Position', [80 80 1300 460]);

subplot(1, 3, 1);
bar(annual_vals, 'FaceColor', [0.85 0.4 0.3]); grid on; hold on;
yline(mean(annual_vals), 'r--', sprintf('均值 %.0f', mean(annual_vals)), 'LineWidth', 1.2);
title(sprintf('全年制氨总量 (吨)\n均=%.0f  最小=%.0f  最大=%.0f', ...
      mean(annual_vals), min(annual_vals), max(annual_vals)));
xlabel('wind+sun 组合'); ylabel('全年制氨量 (吨)');
set(gca, 'XTick', 1:numel(combo_labels_p1), 'XTickLabel', combo_labels_p1, ...
         'XTickLabelRotation', 45);
xlim([0.5, numel(combo_labels_p1)+0.5]);

subplot(1, 3, 2);
bar(cap_vals*100, 'FaceColor', [0.3 0.6 0.85]); grid on; hold on;
yline(mean(cap_vals)*100, 'r--', sprintf('均值 %.2f%%', mean(cap_vals)*100), 'LineWidth', 1.2);
title(sprintf('氢氨年平均产能利用率 (%%)\n均=%.2f  最小=%.2f  最大=%.2f', ...
      mean(cap_vals)*100, min(cap_vals)*100, max(cap_vals)*100));
xlabel('wind+sun 组合'); ylabel('产能利用率 (%)');
set(gca, 'XTick', 1:numel(combo_labels_p1), 'XTickLabel', combo_labels_p1, ...
         'XTickLabelRotation', 45);
xlim([0.5, numel(combo_labels_p1)+0.5]);

subplot(1, 3, 3);
bar(power_vals*100, 'FaceColor', [0.4 0.7 0.4]); grid on; hold on;
yline(mean(power_vals)*100, 'r--', sprintf('均值 %.2f%%', mean(power_vals)*100), 'LineWidth', 1.2);
title(sprintf('风光发电利用率 (%%)\n均=%.2f  最小=%.2f  最大=%.2f', ...
      mean(power_vals)*100, min(power_vals)*100, max(power_vals)*100));
xlabel('wind+sun 组合'); ylabel('发电利用率 (%)');
set(gca, 'XTick', 1:numel(combo_labels_p1), 'XTickLabel', combo_labels_p1, ...
         'XTickLabelRotation', 45);
xlim([0.5, numel(combo_labels_p1)+0.5]);

sgtitle('问题四(1)：离网自治 24 场景产量与利用率指标');


fprintf('\n\n==================== 风光最小装机容量分析 ====================\n');
capacity_results = struct([]);
cap_idx = 0;

for wind_idx = 1:numel(wind_vectors)
    wind = wind_vectors{wind_idx};
    wind_name = wind_names{wind_idx};
    for sun_idx = 1:numel(sun_vectors)
        sun = sun_vectors{sun_idx};
        sun_name = sun_names{sun_idx};

        [Cw_min, Cs_min] = calcMinCapacity(wind, sun, classicload);

        cap_idx = cap_idx + 1;
        capacity_results(cap_idx).wind_name = wind_name;
        capacity_results(cap_idx).sun_name  = sun_name;
        capacity_results(cap_idx).wind_idx  = wind_idx;
        capacity_results(cap_idx).sun_idx   = sun_idx;
        capacity_results(cap_idx).Cw_min    = Cw_min;
        capacity_results(cap_idx).Cs_min    = Cs_min;

        fprintf('%-8s + %-6s：最小风电装机 Cw = %.4g MW，最小光电装机 Cs = %.4g MW\n', ...
                wind_name, sun_name, Cw_min, Cs_min);
    end
end

plotMinCapacity(capacity_results);

function J = objFun(t, wind, sun)
    t = t(:);
  
    wind = wind(:);
    sun = sun(:);

    problem4_cost = 0.75*2*sum(t)*2 ...
                  + 1.5*2*1000*0.2*60000/365/30 ...
                  + 10*2*sum(t)*100 ...
                  + 10*2*sum(t)*150 ...
                  + sum(sun)*64*1000*0.12 ...
                  + sum(wind)*40*1000*0.15;
    J = problem4_cost / (sum(t)*3);
end


function [Cw_min, Cs_min] = calcMinCapacity(wind, sun, classicload)

    wind = wind(:); sun = sun(:); classicload = classicload(:);
    rhs = classicload * 6 + 0.1 * 20.75 * 2;


    f_cap  = [1; 1];
    A_cap  = -[wind, sun];  
    b_cap  = -rhs;
    lb_cap = [0; 0];

    opts_cap = optimoptions('linprog', 'Display', 'off');
    [x_opt, ~, exitflag] = linprog(f_cap, A_cap, b_cap, [], [], lb_cap, [], opts_cap);

    if exitflag == 1
        Cw_min = x_opt(1);
        Cs_min = x_opt(2);
    else
        Cw_min = NaN;
        Cs_min = NaN;
        warning('calcMinCapacity：线性规划无可行解，exitflag=%d。', exitflag);
    end
end

% -------------------------------------------------------------------------
function plotMinCapacity(capacity_results)
    if isempty(capacity_results)
        warning('没有可绘制的装机容量结果。');
        return;
    end

    sort_keys = [[capacity_results.wind_idx]', [capacity_results.sun_idx]'];
    [~, order] = sortrows(sort_keys, [1 2]);
    capacity_results = capacity_results(order);

    combo_labels = arrayfun(@(r) sprintf('%s+%s', r.wind_name, r.sun_name), ...
                            capacity_results, 'UniformOutput', false);
    Cw_vals = [capacity_results.Cw_min];
    Cs_vals = [capacity_results.Cs_min];

    figure('Name', 'problem4：最小风光装机容量', 'Color', 'w');

    subplot(1, 2, 1);
    bar(Cw_vals, 'FaceColor', [0.2 0.6 0.4]);
    grid on;
    title(sprintf('最小风电装机容量 Cw (MW)\n均值=%.4g，最小=%.4g，最大=%.4g', ...
          mean(Cw_vals(~isnan(Cw_vals))), min(Cw_vals), max(Cw_vals)), 'Interpreter', 'none');
    xlabel('wind + sun 组合');
    ylabel('Cw_{min} (MW)', 'Interpreter', 'none');
    set(gca, 'XTick', 1:numel(combo_labels), 'XTickLabel', combo_labels, 'XTickLabelRotation', 45);
    xlim([0.5, numel(combo_labels) + 0.5]);

    subplot(1, 2, 2);
    bar(Cs_vals, 'FaceColor', [0.9 0.6 0.1]);
    grid on;
    title(sprintf('最小光电装机容量 Cs (MW)\n均值=%.4g，最小=%.4g，最大=%.4g', ...
          mean(Cs_vals(~isnan(Cs_vals))), min(Cs_vals), max(Cs_vals)), 'Interpreter', 'none');
    xlabel('wind + sun 组合');
    ylabel('Cs_{min} (MW)', 'Interpreter', 'none');
    set(gca, 'XTick', 1:numel(combo_labels), 'XTickLabel', combo_labels, 'XTickLabelRotation', 45);
    xlim([0.5, numel(combo_labels) + 0.5]);
end

fprintf('\n\n==================== 能源自给与风光利用状况分析 ====================\n');


fprintf('\n---------- 全场景鲁棒最小装机容量 ----------\n');

A_robust  = [];
b_robust  = [];
rhs_const = classicload * 6 + 0.1 * 20.75 * 2;   

for wi = 1:numel(wind_vectors)
    for si = 1:numel(sun_vectors)
        wind_i = wind_vectors{wi};
        sun_i  = sun_vectors{si};
        A_robust = [A_robust; -[wind_i, sun_i]];  
        b_robust = [b_robust; -rhs_const];         
    end
end

f_robust  = [1; 1];
lb_robust = [0; 0];
opts_r    = optimoptions('linprog', 'Display', 'off');
[x_robust, ~, flag_r] = linprog(f_robust, A_robust, b_robust, [], [], lb_robust, [], opts_r);

if flag_r == 1
    Cw_robust = x_robust(1);
    Cs_robust = x_robust(2);
    fprintf('鲁棒最小风电装机 Cw* = %.4f MW\n', Cw_robust);
    fprintf('鲁棒最小光电装机 Cs* = %.4f MW\n', Cs_robust);
    fprintf('合计装机 Cw* + Cs*   = %.4f MW\n', Cw_robust + Cs_robust);
    fprintf('（保证 24 个风光场景的所有时段均不缺电）\n');
else
    Cw_robust = NaN; Cs_robust = NaN;
    warning('鲁棒最小装机求解失败，exitflag=%d', flag_r);
end


figure('Name', 'problem4-1：全场景鲁棒最小装机', 'Color', 'w', ...
       'Position', [200 200 600 450]);
if ~isnan(Cw_robust)
    hb = bar(1, [Cw_robust, Cs_robust], 'stacked'); grid on; hold on;
    hb(1).FaceColor = [0.2 0.6 0.4];
    hb(2).FaceColor = [0.9 0.6 0.1];
    
    % 标注数值
    text(1, Cw_robust/2, sprintf('C_w = %.2f MW', Cw_robust), ...
         'HorizontalAlignment', 'center', 'Color', 'w', 'FontWeight', 'bold');
    text(1, Cw_robust + Cs_robust/2, sprintf('C_s = %.2f MW', Cs_robust), ...
         'HorizontalAlignment', 'center', 'Color', 'w', 'FontWeight', 'bold');
    text(1, Cw_robust + Cs_robust + 0.5, sprintf('合计 = %.2f MW', Cw_robust + Cs_robust), ...
         'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    
    % 与题目给定装机对比
    yline(40, 'r--', '题目风电 40 MW', 'LineWidth', 1.2);
    yline(40 + 64, 'b--', '题目合计 104 MW', 'LineWidth', 1.2);
    
    set(gca, 'XTick', 1, 'XTickLabel', {'全场景鲁棒最小装机'});
    ylabel('装机容量 (MW)');
    legend({'C_w (风电)', 'C_s (光电)'}, 'Location', 'best');
    title('能源自治最小风光装机（保证 24 场景全可行）');
end




set(0, 'DefaultAxesFontName', 'SimHei');
set(0, 'DefaultTextFontName', 'SimHei');
set(0, 'DefaultAxesFontSize', 10);

P_w  = 40;     P_s = 64;     P_L = 6;
P_h2 = 20.75 * 2;            
n_t  = 24;                    


fprintf('\n==================== 问题四(1) 自给状况详细分析 ====================\n');
fprintf('%-15s %-10s %-10s %-10s %-10s %-10s %-10s\n', ...
        '场景', '日发电', '日总用电', '风光利用率', '自给率', '弃电率', '缺电时段');

p41_extra = struct([]);
for k = 1:numel(problem4_results)
    r = problem4_results(k);
    wind = wind_vectors{r.wind_idx};
    sun  = sun_vectors{r.sun_idx};
    
    gen_hour    = wind*P_w + sun*P_s;            % 各时段发电 (MW)
    load_h2_hr  = r.t_opt * P_h2;                % 各时段制氢氨耗电
    load_other  = classicload * P_L;             % 各时段常规负荷
    load_total  = load_h2_hr + load_other;       % 各时段总用电
    

    surplus     = gen_hour - load_total;
    curt_hour   = max(surplus, 0);               % 弃风弃光
    deficit_hr  = max(-surplus, 0);              % 理论上调度可行时应≈0
    
    daygen      = sum(gen_hour);
    dayload     = sum(load_total);
    daycurt     = sum(curt_hour);
    daydef      = sum(deficit_hr);
    
    selfsuff    = (dayload - daydef) / dayload;          % 自给率
    power_use   = (daygen - daycurt) / daygen;           % 风光利用率
    curt_rate   = daycurt / daygen;                       % 弃电率
    
    p41_extra(k).case_idx   = r.case_idx;
    p41_extra(k).wind_name  = r.wind_name;
    p41_extra(k).sun_name   = r.sun_name;
    p41_extra(k).wind_idx   = r.wind_idx;
    p41_extra(k).sun_idx    = r.sun_idx;
    p41_extra(k).daygen     = daygen;
    p41_extra(k).dayload    = dayload;
    p41_extra(k).daycurt    = daycurt;
    p41_extra(k).curt_rate  = curt_rate;
    p41_extra(k).self_suff  = selfsuff;
    p41_extra(k).power_use  = power_use;
    p41_extra(k).gen_hour   = gen_hour;
    p41_extra(k).load_total = load_total;
    p41_extra(k).curt_hour  = curt_hour;
    p41_extra(k).low_hours  = sum(r.t_opt < 0.15);       % 接近最低负荷的时段数
    
    fprintf('%-15s %-10.1f %-10.1f %-10.2f%% %-10.2f%% %-10.2f%% %-10d\n', ...
            sprintf('%s+%s', r.wind_name, r.sun_name), ...
            daygen, dayload, power_use*100, selfsuff*100, curt_rate*100, ...
            p41_extra(k).low_hours);
end

% ---------- 排序后绘 4 联图 ----------
sort_keys = [[p41_extra.wind_idx]', [p41_extra.sun_idx]'];
[~, ord]  = sortrows(sort_keys, [1 2]);
p41_extra = p41_extra(ord);
labels = arrayfun(@(r) sprintf('%s+%s', r.wind_name, r.sun_name), p41_extra, 'UniformOutput', false);

vals_self = [p41_extra.self_suff] * 100;
vals_use  = [p41_extra.power_use] * 100;
vals_curt = [p41_extra.curt_rate] * 100;
vals_gen  = [p41_extra.daygen];
vals_load = [p41_extra.dayload];

figure('Name', '问题四(1)：能源自给与风光利用分析', 'Color', 'w', ...
       'Position', [60 60 1400 700]);


subplot(2,2,1);
bar_data = [vals_self; vals_use]';
b_h = bar(bar_data, 'grouped'); grid on; hold on;
b_h(1).FaceColor = [0.20 0.60 0.40];   % 自给率
b_h(2).FaceColor = [0.30 0.55 0.85];   % 利用率
yline(100, 'k--', '100%', 'LineWidth', 1.0);
yline(mean(vals_self), ':', sprintf('自给均值%.1f%%', mean(vals_self)), ...
      'Color', [0.20 0.60 0.40], 'LineWidth', 1.2);
title('能源自给率 vs 风光发电利用率');
ylabel('百分比 (%)'); xlabel('风光场景');
legend({'能源自给率', '风光利用率'}, 'Location', 'southwest', 'FontSize', 9);
set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels, 'XTickLabelRotation', 45);
xlim([0.5, numel(labels)+0.5]); ylim([0 110]);


subplot(2,2,2);
bar(vals_curt, 'FaceColor', [0.90 0.40 0.20]); grid on; hold on;
yline(mean(vals_curt), 'r--', sprintf('均值 %.2f%%', mean(vals_curt)), 'LineWidth', 1.2);
title(sprintf('风光弃电率（均=%.2f%%，最大=%.2f%%）', mean(vals_curt), max(vals_curt)));
xlabel('风光场景'); ylabel('弃电率 (%)');
set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels, 'XTickLabelRotation', 45);
xlim([0.5, numel(labels)+0.5]);


subplot(2,2,3);
plot(vals_gen, '-o', 'LineWidth', 1.5, 'Color', [0.30 0.55 0.85], ...
     'MarkerFaceColor', [0.30 0.55 0.85], 'DisplayName', '日发电量');
hold on; grid on;
plot(vals_load, '-s', 'LineWidth', 1.5, 'Color', [0.85 0.40 0.20], ...
     'MarkerFaceColor', [0.85 0.40 0.20], 'DisplayName', '日用电量');

fill_x = [1:numel(vals_gen), numel(vals_gen):-1:1];
fill_y = [vals_gen, fliplr(vals_load)];
fill(fill_x, fill_y, [0.7 0.85 0.7], 'FaceAlpha', 0.3, ...
     'EdgeColor', 'none', 'DisplayName', '可弃电缓冲');
title('日发电量 vs 日用电量（绿色填充=可消纳冗余）');
xlabel('风光场景'); ylabel('电量 (MW·h)');
legend('Location', 'best', 'FontSize', 9);
set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels, 'XTickLabelRotation', 45);
xlim([0.5, numel(labels)+0.5]);


subplot(2,2,4);
margin = (vals_gen - vals_load) ./ vals_load * 100;   % 发电相对用电的盈余 (%)
scatter(margin, vals_curt, 80, [P_w*0 0.30 0.70], 'filled'); grid on; hold on;
for k = 1:numel(labels)
    text(margin(k)+0.5, vals_curt(k), labels{k}, 'FontSize', 7);
end
xlabel('日发电盈余率 (发电-用电)/用电 (%)'); ylabel('弃电率 (%)');
title('发电盈余 vs 弃电率（接近线性正相关）');

p = polyfit(margin, vals_curt, 1);
xx = linspace(min(margin), max(margin), 50);
plot(xx, polyval(p, xx), 'r--', 'LineWidth', 1.2, ...
     'DisplayName', sprintf('y = %.3fx + %.2f', p(1), p(2)));
legend('show', 'Location', 'northwest');

sgtitle('问题四(1)：离网无储能场景下园区能源自给与风光利用分析');
saveas(gcf, 'problem4_1_self_sufficiency.png');



