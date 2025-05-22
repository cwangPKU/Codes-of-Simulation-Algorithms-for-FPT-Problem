function plot_brownian_density(W1, W2, in_region, a1, a2, b1, b2)
    % 提取最终位置
    final_W1 = W1(:, end);
    final_W2 = W2(:, end);
    
    % 仅保留未越界路径的最终位置
    valid_final = [final_W1(in_region), final_W2(in_region)];
    survival_prob = mean(in_region);
    
    % 设置网格
    grid_size = 50;
    x_edges = linspace(a1, b1, grid_size);
    y_edges = linspace(a2, b2, grid_size);
    
    % 计算2D直方图（绝对计数）
    [counts, x_centers, y_centers] = histcounts2(...
        valid_final(:,1), valid_final(:,2), x_edges, y_edges);
    
    % 关键修正：密度归一化
    bin_area = (x_edges(2)-x_edges(1)) * (y_edges(2)-y_edges(1)); % 每个bin的面积
    density = counts / (bin_area * size(W1,1)); % 密度=计数/(总面积*总路径数)
    
    % 验证积分=留存概率
    fprintf('Survival Probability: %.4f\n', survival_prob);
    fprintf('Intergration of Density: %.4f\n', sum(density(:)) * bin_area);
    
    % 绘图
    figure;
    imagesc(x_centers, y_centers, density');
    axis xy; colorbar; hold on;
    rectangle('Position', [a1,a2,b1-a1,b2-a2], 'EdgeColor', 'r', 'LineWidth', 1.5);
    title(sprintf('Monte Carlo Method (Survival Prob=%.2f%%)', survival_prob*100));
    xlabel('W1(T)'); ylabel('W2(T)');
    colormap('jet');
end