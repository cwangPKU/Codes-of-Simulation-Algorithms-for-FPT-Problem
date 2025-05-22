function [prob_a, prob_b, prob_in, W_T] = bm1d_sim(a, b, T, K, N)

    dt = T / N;     % 步长
    
    % 生成 K 条布朗运动路径                     % 固定随机种子（可选）
    dW = sqrt(dt) * randn(K, N);        % K×N 的随机增量矩阵
    W = [zeros(K,1), cumsum(dW, 2)];    % 沿行累积求和，初始位置为0
    
    hit_a = (W <= a);
    hit_b = (W >= b);
    
    [~, idx_a] = max(hit_a, [], 2); 
    [~, idx_b] = max(hit_b, [], 2);
    
    idx_a(idx_a==1) = Inf;
    idx_b(idx_b==1) = Inf;
    
    [min_id, first_touch] = min([idx_a, idx_b], [], 2);
    
    count_a = sum ((min_id < Inf) & (first_touch == 1));
    count_b = sum ((min_id < Inf) & (first_touch == 2));
    count_in = sum(min_id == Inf);
    
    prob_a = count_a / K;
    prob_b = count_b / K;
    prob_in = count_in / K;
    
    W_T = W(min_id == Inf, N+1);
end