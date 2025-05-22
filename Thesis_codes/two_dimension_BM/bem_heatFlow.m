function [heatFlow, points_x, points_y, S] = bem_heatFlow(a1, a2, b1, b2, rho, T, m, n)

    points_x = [];
    points_y = [];
    S = [];
    
    points1 = linspace(a1, b1, m+2);
    points1 = points1(2: end-1);
    points2 = linspace(a2, b2, m+2);
    points2 = points2(2: end-1);
    
    points_x = [points_x, points1];
    points_y = [points_y, (a2 - rho * points1) ./ sqrt(1 - rho^2)];
    S = [S, (b1- a1) / (m * sqrt(1-rho^2)) * ones(1, m)];
    
    points_x = [points_x, points1];
    points_y = [points_y, (b2 - rho * points1) ./ sqrt(1 - rho^2)];
    S = [S, (b1- a1) / (m * sqrt(1-rho^2)) * ones(1, m)];
    
    points_x = [points_x, a1 * ones(1, m)];
    points_y = [points_y, (points2 - rho * a1) ./ sqrt(1 - rho^2)];
    S = [S, (b2- a2) / (m * sqrt(1-rho^2)) * ones(1, m)];
    
    points_x = [points_x, b1 * ones(1, m)];
    points_y = [points_y, (points2 - rho * b1) ./ sqrt(1 - rho^2)];
    S = [S, (b2- a2) / (m * sqrt(1-rho^2)) * ones(1, m)];
    

    x_mat = points_x' - points_x;
    y_mat = points_y' - points_y;

    dt = T / n;

    heatFlow = zeros(n, m*4);
    A = (dt / 2) .* S .* G(x_mat, y_mat, dt);


    for k = 1:n
        b = zeros(4*m, 1);
        for a = 0:k-2
            b = b + sum(G(x_mat, y_mat, (k-a)*dt) .* heatFlow(a+1, :) .* S, 2);
        end
        b = -1 * b * dt / 2;
        b = b - G(points_x', points_y', dt * k);
        
        heatFlow(k, :) = (A \ b)';

    end
end