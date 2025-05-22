function [W1, W2, in_region, exit_time, exit_pos] = bm2d_sim(a1, a2, b1, b2, rho, T, n, N)
    % 参数校验（同上，此处省略）
    
    % 生成布朗运动路径（同上）
    dt = T/n;
    t = linspace(0, T, n+1);
    Z1 = randn(N, n);
    Z2 = randn(N, n);
    dW1 = sqrt(dt) * Z1;
    dW2 = rho*sqrt(dt)*Z1 + sqrt(1-rho^2)*sqrt(dt)*Z2;
    W1 = [zeros(N,1), cumsum(dW1,2)];
    W2 = [zeros(N,1), cumsum(dW2,2)];

    % 初始化输出变量
    [in_region, exit_time, exit_pos] = deal(true(N,1), nan(N,1), nan(N,2));
    out_of_bounds = (W1 <= a1) | (W1 >= b1) | (W2 <= a2) | (W2 >= b2);

    for k = 1:N
        % 找到首次越界的时间步
        first_out = find(out_of_bounds(k,2:end), 1); % 从第2步开始检查
        if ~isempty(first_out)
            in_region(k) = false;
            exit_time(k) = t(first_out+1); % 越界时间
            
            % 越界前的最后一个内部点 (P1)
            P1 = [W1(k,first_out), W2(k,first_out)];
            % 越界后的外部点 (P2)
            P2 = [W1(k,first_out+1), W2(k,first_out+1)];
            
            % 计算连线与边界的交点
            exit_pos(k,:) = find_boundary_intersection(P1, P2, a1, a2, b1, b2);
        end
    end
end

% 辅助函数：计算线段与矩形边界的交点
function intersection = find_boundary_intersection(P1, P2, a1, a2, b1, b2)
    % 提取坐标
    x1 = P1(1); y1 = P1(2);
    x2 = P2(1); y2 = P2(2);
    
    % 检查与四条边界的交点（按顺序：左、右、下、上）
    boundaries = {[a1, a2, a1, b2], [b1, a2, b1, b2], [a1, a2, b1, a2], [a1, b2, b1, b2]};
    
    for i = 1:4
        b = boundaries{i};
        % 线段参数方程：P = P1 + t*(P2-P1), t∈[0,1]
        % 边界线参数：x3,y3 -> x4,y4
        x3 = b(1); y3 = b(2);
        x4 = b(3); y4 = b(4);
        
        % 计算两线段的交点
        denominator = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4);
        if abs(denominator) < 1e-10
            continue; % 平行或共线，跳过
        end
        
        t = ((x1-x3)*(y3-y4) - (y1-y3)*(x3-x4)) / denominator;
        u = -((x1-x2)*(y1-y3) - (y1-y2)*(x1-x3)) / denominator;
        
        if t >= 0 && t <= 1 && u >= 0 && u <= 1
            intersection = [x1 + t*(x2-x1), y1 + t*(y2-y1)];
            return;
        end
    end
    
    % 理论上应该总能找到交点，此处为容错处理
    intersection = P2;
end