function f = laplace1d(t, a, b, N, x)
    
    [prob_a_theo, prob_b_theo, prob_in_theo] = prob1d(t, a, b ,20);
    F = @(s) exp(0.5 * s.^2 * t) - exp(a * s) * prob_a_theo - exp(b * s) * prob_b_theo;
    f = F(x);
    %f = stehfest_vec(F, x, N);
    
end

function f = talbot_vec(F, t, N)
    % F: Laplace 变换函数句柄（需支持矩阵输入，例如 F = @(s) 1./(s.^2 + 1)）
    % t: 时间点（标量或向量）
    % N: 积分节点数
    t = t(:); % 转换为列向量 (M x 1)
    theta = linspace(-pi, pi, N); % 生成积分节点 (1 x N)
    sigma = 0.2 ./ t; % 路径参数 (M x 1)
    
    % 构建复数积分路径 z (M x N)
    z = sigma .* t + 1i * theta; 
    s = z ./ t; % Laplace 变量 s (M x N)
    
    % 计算积分项
    dsdt = (1i * sigma - 1i * z) ./ t; % 导数项 (M x N)
    integrand = F(s) .* exp(z) .* dsdt; % 被积函数 (M x N)
    
    % 求和并计算最终结果
    sum_vals = sum(integrand, 2); % 按行求和 (M x 1)
    d_theta = theta(2) - theta(1); % 积分步长
    f = real(sum_vals * d_theta / (2*pi*1i)); 
    f = reshape(f, size(t)); % 保持输出形状与输入 t 一致
end


function f = stehfest_vec(F, t, N)
    % F: 支持向量化输入的 Laplace 变换函数句柄（例如 F = @(s) 1./(s+1)）
    % t: 时间点（可以是标量或向量）
    % N: 权重项数（建议偶数，如 6-18）
    V = stehfest_weights(N);       % 计算权重系数
    ln2 = log(2);
    k = 1:N;                        % 生成 k 向量 [1, 2, ..., N]
    s = k' * ln2 ./ t(:)';          % 构建 s 矩阵（维度：length(t) × N）
    F_values = F(s);                % 向量化计算 F(s)
    sum_vals = sum(F_values .* V(:), 1);     % 矩阵乘法求和（避免显式循环）
    f = (ln2 ./ t) .* sum_vals;      % 最终结果（支持 t 为向量）
end

% 保持权重计算函数不变（因 N 较小，向量化收益有限）
function V = stehfest_weights(N)
    V = zeros(1, N);
    for k = 1:N
        j_range = floor((k+1)/2):min(k, N/2);
        j = j_range(:);
        terms = (j.^(N/2)) .* factorial(2*j) ./ ...
            (factorial(N/2 - j) .* factorial(j) .* factorial(j-1) .* ...
             factorial(k - j) .* factorial(2*j - k));
        sum_val = sum(terms);
        V(k) = (-1)^(k + N/2) * sum_val;
    end
end