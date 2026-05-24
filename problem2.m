n = length(x);
opts = optimoptions('ga', 'Display', 'iter');

beq_values = [24, 21, 18, 15, 12];
Aeq = ones(1, n);

problem2_results = struct([]);

for i = 1:numel(beq_values)
    beq = beq_values(i);

    fprintf('\n========== beq = %g, sum(t) = %g ==========\n', beq, beq);

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
    problem2_sell = problem2_buy;
    problem2_buy(problem2_buy<0) = 0;
    problem2_sell(problem2_sell>0) = 0;
    problem2_buycost = sum(problem2_buy .* fenshijijia);
    problem2_dayload = sum(problem2_load);
    problem2_daygenerate = sum(problem2_generate);
    problem2_daybuy = sum(problem2_buy);
    problem2_daysell = abs(sum(problem2_sell));
    problem2_xnyzf = (problem2_dayload-problem2_daysell-problem2_daybuy)/problem2_daygenerate;
    problem2_zydl = (problem2_daygenerate-problem2_daysell)/problem2_dayload;
    problem2_xnysw = problem2_daysell/problem2_daygenerate;
    fprintf('新能源自发自用电量占总可用发电量比例 = %.6g\n', problem2_xnyzf);
    fprintf('总用电量中绿电占比 = %.6g\n', problem2_zydl);
    fprintf('新能源上网电量比例 = %.6g\n', problem2_xnysw);

    problem2_results(i).beq = beq;
    problem2_results(i).t_opt = t_opt;
    problem2_results(i).J_min = J_min;
    problem2_results(i).nh3load = problem2_nh3load;
    problem2_results(i).load = problem2_load;
    problem2_results(i).generate = problem2_generate;
    problem2_results(i).buy = problem2_buy;
    problem2_results(i).sell = problem2_sell;
    problem2_results(i).buycost = problem2_buycost;
    problem2_results(i).dayload = problem2_dayload;
    problem2_results(i).daygenerate = problem2_daygenerate;
    problem2_results(i).daybuy = problem2_daybuy;
    problem2_results(i).daysell = problem2_daysell;
    problem2_results(i).xnyzf = problem2_xnyzf;
    problem2_results(i).zydl = problem2_zydl;
    problem2_results(i).xnysw = problem2_xnysw;
end

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

wind_vectors = {wind1(:), wind2(:), wind3(:), wind4(:), wind5(:), wind6(:)};
wind_names = {'wind1', 'wind2', 'wind3', 'wind4', 'wind5', 'wind6'};
sun_vectors = {sun1(:), sun2(:), sun3(:), sun4(:)};
sun_names = {'sun1', 'sun2', 'sun3', 'sun4'};

classicload = classicload(:);
fenshijijia = fenshijijia(:);
n = numel(wind_vectors{1});


opts = optimoptions('ga', 'Display', 'iter');


beq_values = [24, 21, 18, 15, 12];
Aeq = ones(1, n);

problem2_results = struct([]);
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

            t_opt = ga(@(t) objFun(t, classicload, wind, sun, fenshijijia), ...
                       n, [], [], Aeq, beq, zeros(n,1), ones(n,1), [], 1:n, opts);
            t_opt = t_opt(:);
            J_min = objFun(t_opt, classicload, wind, sun, fenshijijia);
            fprintf('最小目标值 J = %.6g\n', J_min);
            disp('最优 t 向量：'); disp(t_opt);

            metrics = calcProblem2Metrics(t_opt, classicload, wind, sun, fenshijijia);
                                         

            problem2_nh3load = metrics.nh3load;
            problem2_load = metrics.load;
            problem2_generate = metrics.generate;
            problem2_buy = metrics.buy;
            problem2_sell = metrics.sell;
            problem2_buycost = metrics.buycost;
            problem2_dayload = metrics.dayload;
            problem2_daygenerate = metrics.daygenerate;
            problem2_daybuy = metrics.daybuy;
            problem2_daysell = metrics.daysell;
            problem2_xnyzf = metrics.xnyzf;
            problem2_zydl = metrics.zydl;
            problem2_xnysw = metrics.xnysw;

            fprintf('新能源自发自用电量占总可用发电量比例 = %.6g\n', problem2_xnyzf);
            fprintf('总用电量中绿电占比 = %.6g\n', problem2_zydl);
            fprintf('新能源上网电量比例 = %.6g\n', problem2_xnysw);

            result_idx = result_idx + 1;
            problem2_results(result_idx).case_idx = result_idx;
            problem2_results(result_idx).wind_idx = wind_idx;
            problem2_results(result_idx).wind_name = wind_name;
            problem2_results(result_idx).sun_idx = sun_idx;
            problem2_results(result_idx).sun_name = sun_name;
            problem2_results(result_idx).beq = beq;
            problem2_results(result_idx).t_opt = t_opt;
            problem2_results(result_idx).J_min = J_min;
            problem2_results(result_idx).nh3load = problem2_nh3load;
            problem2_results(result_idx).load = problem2_load;
            problem2_results(result_idx).generate = problem2_generate;
            problem2_results(result_idx).buy = problem2_buy;
            problem2_results(result_idx).sell = problem2_sell;
            problem2_results(result_idx).buycost = problem2_buycost;
            problem2_results(result_idx).dayload = problem2_dayload;
            problem2_results(result_idx).daygenerate = problem2_daygenerate;
            problem2_results(result_idx).daybuy = problem2_daybuy;
            problem2_results(result_idx).daysell = problem2_daysell;
            problem2_results(result_idx).xnyzf = problem2_xnyzf;
            problem2_results(result_idx).zydl = problem2_zydl;
            problem2_results(result_idx).xnysw = problem2_xnysw;
        end
    end
end

plotProblem2Distributions(problem2_results, beq_values);

function metrics = calcProblem2Metrics(t, classicload, wind, sun, fenshijijia )
                                     
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


function plotProblem2Distributions(problem2_results, beq_values)
    metric_fields = {'xnyzf', 'zydl', 'xnysw', 'J_min', 'daybuy', 'daysell'};
    metric_names = {'problem2\_xnyzf', 'problem2\_zydl', 'problem2\_xnysw', ...
                    'J', 'problem2\_daybuy', 'problem2\_daysell'};

    for beq_idx = 1:numel(beq_values)
        beq = beq_values(beq_idx);
        beq_mask = [problem2_results.beq] == beq;
        beq_results = problem2_results(beq_mask);

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



assert(exist('problem2_results','var')==1, 'problem2_results 不存在');


n_records = numel(problem2_results);
if n_records ~= 120
    error(['problem2_results 当前有 %d 条记录，期望 120 条（24场景×5产量）。\n' ...
           '请先运行 24 种 wind+sun 组合的双层循环优化代码。'], n_records);
end


fns = fieldnames(problem2_results);
if ~ismember('wind_idx', fns)
    for k = 1:n_records
        problem2_results(k).wind_idx = floor((k-1)/(4*5)) + 1;
        problem2_results(k).sun_idx  = mod(floor((k-1)/5), 4) + 1;
        problem2_results(k).wind_name = sprintf('wind%d', problem2_results(k).wind_idx);
        problem2_results(k).sun_name  = sprintf('sun%d',  problem2_results(k).sun_idx);
    end
    warning('problem2_results 中无 wind_idx 字段，已按 case_idx 顺序自动补齐。请核对场景顺序是否为 wind 外层、sun 内层。');
end

all_wind_idx = [problem2_results.wind_idx];
all_sun_idx  = [problem2_results.sun_idx];

n_wind = max(all_wind_idx);   % 6
n_sun  = max(all_sun_idx);    % 4
n_scen = n_wind * n_sun;      % 24
days_per_scen = 15;
n_days_year   = n_scen * days_per_scen;   % 360

best_J        = zeros(n_scen,1);
best_beq      = zeros(n_scen,1);
best_xnyzf    = zeros(n_scen,1);
best_zydl     = zeros(n_scen,1);
best_xnysw    = zeros(n_scen,1);
best_label    = strings(n_scen,1);
best_prod_ton = zeros(n_scen,1);

beq_to_ton = @(b) b * 1.5 * 2;   % 日产氨吨数 = sum(t)*1.5*2

scen_id = 0;
for wi = 1:n_wind
    for si = 1:n_sun
        scen_id = scen_id + 1;
        mask = (all_wind_idx==wi) & (all_sun_idx==si);
        sub = problem2_results(mask);
        if isempty(sub)
            error('场景 wind%d+sun%d 在 problem2_results 中找不到记录', wi, si);
        end
        [J_min_val, k] = min([sub.J_min]);
        best_J(scen_id)     = J_min_val;
        best_beq(scen_id)   = sub(k).beq;
        best_xnyzf(scen_id) = sub(k).xnyzf;
        best_zydl(scen_id)  = sub(k).zydl;
        best_xnysw(scen_id) = sub(k).xnysw;
        best_label(scen_id) = sprintf('w%d+s%d', wi, si);
        best_prod_ton(scen_id) = beq_to_ton(sub(k).beq);
    end
end

% 绿电指标合格性
ok_xnyzf = best_xnyzf > 0.60;
ok_zydl  = best_zydl  > 0.30;
ok_xnysw = best_xnysw < 0.20;
ok_all   = ok_xnyzf & ok_zydl & ok_xnysw;
n_pass_each = [ok_xnyzf, ok_zydl, ok_xnysw];
n_pass      = sum(n_pass_each, 2);   % 每个场景通过几项（0/1/2/3）
ok_part     = (n_pass>=1) & (n_pass<=2);
ok_none     = (n_pass==0);

n_full = sum(ok_all);  n_part = sum(ok_part);  n_none = sum(ok_none);

year_J     = repelem(best_J,        days_per_scen);
year_prod  = repelem(best_prod_ton, days_per_scen);

total_year_cost = sum(year_J .* year_prod);
total_year_ton  = sum(year_prod);
avg_year_cost   = total_year_cost / total_year_ton;

fig = figure('Name','问题二(2) 全年吨氨成本分布','Color','w',...
             'Position',[100 100 1280 780]);
set(fig,'DefaultAxesFontName','SimHei','DefaultTextFontName','SimHei');

c_full = [0.20 0.65 0.35];
c_part = [0.95 0.70 0.20];
c_none = [0.85 0.30 0.30];
c_line = [0.20 0.30 0.55];


ax1 = subplot(2,2,[1 2]); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');

x = 1:n_days_year;
y_lo = min(year_J)*0.95;  y_hi = max(year_J)*1.05;

for s = 1:n_scen
    x_start = (s-1)*days_per_scen + 0.5;
    x_end   =  s   *days_per_scen + 0.5;
    if ok_all(s),       fc = c_full;
    elseif ok_part(s),  fc = c_part;
    else,               fc = c_none;
    end
    patch(ax1,[x_start x_end x_end x_start], ...
          [y_lo y_lo y_hi y_hi], ...
          fc,'FaceAlpha',0.18,'EdgeColor','none','HandleVisibility','off');
end

plot(ax1, x, year_J, '-', 'Color', c_line, 'LineWidth', 1.6);
scatter(ax1, x, year_J, 14, c_line, 'filled', 'MarkerFaceAlpha', 0.6, 'HandleVisibility','off');

yline(ax1, avg_year_cost, '--', sprintf('全年加权平均 = %.2f 元/吨', avg_year_cost), ...
      'Color',[0.4 0.4 0.4], 'LineWidth', 1.3, ...
      'LabelHorizontalAlignment','left','FontSize', 10);

for s = 1:n_scen-1
    xline(ax1, s*days_per_scen+0.5, ':', 'Color',[0.75 0.75 0.75], 'HandleVisibility','off');
end


for s = 1:n_scen
    text(ax1, (s-0.5)*days_per_scen+0.5, y_hi*0.995, best_label(s), ...
        'HorizontalAlignment','center','VerticalAlignment','top', ...
        'FontSize',7,'Color',[0.3 0.3 0.3]);
end

xlabel(ax1,'全年天数（每场景代表 15 天，共 360 天）','FontSize',11);
ylabel(ax1,'吨氨成本 J （元/吨）','FontSize',11);
title(ax1, sprintf('园区全年吨氨成本分布曲线  |  全年总吨氨成本 = %.4e 元  |  全年总产量 = %.1f 吨', ...
       total_year_cost, total_year_ton), 'FontSize',12);
xlim(ax1,[0.5 n_days_year+0.5]);
ylim(ax1,[y_lo y_hi]);

h1 = patch(nan,nan,c_full,'FaceAlpha',0.25,'EdgeColor','none');
h2 = patch(nan,nan,c_part,'FaceAlpha',0.25,'EdgeColor','none');
h3 = patch(nan,nan,c_none,'FaceAlpha',0.25,'EdgeColor','none');
h4 = plot(ax1,nan,nan,'-','Color',c_line,'LineWidth',1.6);
legend(ax1,[h4 h1 h2 h3], ...
       {'吨氨成本','三项指标全满足','部分指标满足','三项指标全不满足'}, ...
       'Location','best','FontSize',10);


ax2 = subplot(2,2,3); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
bar_colors = repmat(c_part,n_scen,1);
bar_colors(ok_all,:)  = repmat(c_full, n_full, 1);
bar_colors(ok_none,:) = repmat(c_none, n_none, 1);

b = bar(ax2, 1:n_scen, best_J, 'FaceColor','flat','EdgeColor',[0.3 0.3 0.3]);
b.CData = bar_colors;

[J_lo, i_lo] = min(best_J);
[J_hi, i_hi] = max(best_J);
text(ax2, i_lo, J_lo, sprintf('↓最低\n%.1f', J_lo), ...
     'HorizontalAlignment','center','VerticalAlignment','bottom',...
     'FontSize',9,'Color',[0.1 0.5 0.2],'FontWeight','bold');
text(ax2, i_hi, J_hi, sprintf('↑最高\n%.1f', J_hi), ...
     'HorizontalAlignment','center','VerticalAlignment','bottom',...
     'FontSize',9,'Color',[0.7 0.2 0.2],'FontWeight','bold');

set(ax2,'XTick',1:n_scen,'XTickLabel',cellstr(best_label),...
    'XTickLabelRotation',55,'FontSize',8);
xlabel(ax2,'风光组合场景','FontSize',10);
ylabel(ax2,'最优吨氨成本（元/吨）','FontSize',10);
title(ax2,'24种风光场景下最优吨氨成本','FontSize',11);
xlim(ax2,[0.5 n_scen+0.5]);


ax3 = subplot(2,2,4);
days_full = n_full * days_per_scen;
days_part = n_part * days_per_scen;
days_none = n_none * days_per_scen;
pie_data   = [days_full, days_part, days_none];
pie_labels = {sprintf('全满足\n%d天 (%.1f%%)', days_full, days_full/n_days_year*100), ...
              sprintf('部分满足\n%d天 (%.1f%%)', days_part, days_part/n_days_year*100), ...
              sprintf('全不满足\n%d天 (%.1f%%)', days_none, days_none/n_days_year*100)};
valid = pie_data > 0;
colors_pie = [c_full; c_part; c_none];
colors_pie = colors_pie(valid,:);
ph = pie(ax3, pie_data(valid), pie_labels(valid));
patches_idx = 1;
for kk = 1:numel(ph)
    if isa(ph(kk),'matlab.graphics.primitive.Patch')
        ph(kk).FaceColor = colors_pie(patches_idx,:);
        ph(kk).EdgeColor = 'w';
        ph(kk).LineWidth = 1.2;
        patches_idx = patches_idx + 1;
    elseif isa(ph(kk),'matlab.graphics.primitive.Text')
        ph(kk).FontSize = 9.5;
        ph(kk).FontName = 'SimHei';
    end
end
title(ax3,'全年绿电直连指标合格情况（按天数）','FontSize',11);

sgtitle('问题二（2）：基于离散制氨调节的全年运行结果分析','FontSize',13,'FontWeight','bold');

fprintf('\n========== 问题二(2) 全年统计汇总 ==========\n');
fprintf('全年天数         : %d 天（24场景 × 15天）\n', n_days_year);
fprintf('全年总产氨量     : %.2f 吨\n', total_year_ton);
fprintf('全年总吨氨成本   : %.4e 元\n', total_year_cost);
fprintf('全年加权平均成本 : %.2f 元/吨\n', avg_year_cost);
fprintf('最低吨氨成本场景 : %s ，J = %.2f 元/吨\n', best_label(i_lo), J_lo);
fprintf('最高吨氨成本场景 : %s ，J = %.2f 元/吨\n', best_label(i_hi), J_hi);
fprintf('绿电指标 全满足  : %d 场景 / %d 天 (%.1f%%)\n', n_full, days_full, days_full/n_days_year*100);
fprintf('绿电指标 部分满足: %d 场景 / %d 天 (%.1f%%)\n', n_part, days_part, days_part/n_days_year*100);
fprintf('绿电指标 全不满足: %d 场景 / %d 天 (%.1f%%)\n', n_none, days_none, days_none/n_days_year*100);

T = table((1:n_scen)', best_label, best_prod_ton, best_beq, best_J, ...
          best_xnyzf, best_zydl, best_xnysw, ok_all, ...
    'VariableNames', {'场景','风光组合','日产氨_吨','sum_t','吨氨成本',...
                      '自发自用比','绿电占比','上网比','三项全合格'});
disp(T);
