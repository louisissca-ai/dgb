% 横坐标：0 到 23 的整数
x = (0:23)';

% 曲线 1 的纵坐标（填入 24 个数值）
y1 = [
];

% 曲线 2 的纵坐标（填入 24 个数值）
y2 = [
];

% 曲线 3 的纵坐标（填入 24 个数值）
y3 = [
];

% 曲线 4 的纵坐标（填入 24 个数值）
y4 = [
];

% 绘制线型图
figure;
hold on;
plot(x, y1, '-o', 'DisplayName', '曲线1');
plot(x, y2, '-s', 'DisplayName', '曲线2');
plot(x, y3, '-^', 'DisplayName', '曲线3');
plot(x, y4, '-d', 'DisplayName', '曲线4');
hold off;

xlabel('小时');
ylabel('数值');
title('线型图');

xticks(0:23);
legend('show');
grid on;
