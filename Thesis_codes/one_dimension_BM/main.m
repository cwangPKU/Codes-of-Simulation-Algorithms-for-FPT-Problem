a = -3;
b = 2;
t = 2;

K = 1e5;
N = 1e4;

[prob_a, prob_b, prob_in, W_T] = bm1d_sim(a, b, t, K, N);

[prob_a_theo, prob_b_theo, prob_in_theo] = prob1d(t, a, b ,N);

x = linspace(a, b, 1e5);
ans_reflection = reflection1d(x, t, a, b, -1e3, 1e3);
ans_green = green1d(x, t, a, b, 1e3);


[counts, edges] = histcounts(W_T);
binWidth = edges(2) - edges(1); % 分箱宽度（假设等宽）
originalArea = sum(counts) * binWidth;
scalingFactor = prob_in / originalArea;
scaledCounts = counts * scalingFactor;

% 绘制直方图
figure;
bar(edges(1:end-1), scaledCounts, 1, 'EdgeColor', 'none', 'DisplayName', 'Discrete Simulation');
hold on;
plot_reflection = plot(x, ans_reflection, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Reflection Principle');
plot_green = plot(x, ans_green, 'g--', 'LineWidth', 2.5, 'DisplayName', 'PDE Method');

% 图形装饰
xlabel('Position (x)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Conditional Probability Density', 'FontSize', 12, 'FontWeight', 'bold');
title('Conditional Probability Density', 'FontSize', 14);
legend([plot_reflection, plot_green], 'Location', 'best'); % 自动选择最佳图例位置
grid on;
box on;

disp([prob_a, prob_b, prob_in]);
disp([prob_a_theo, prob_b_theo, prob_in_theo]);