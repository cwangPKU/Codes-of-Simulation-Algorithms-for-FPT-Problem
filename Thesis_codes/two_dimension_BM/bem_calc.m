function density = bem_calc(x, y, rho, T, heatFlow, points_x, points_y, S)

x_new = x;
y_new = (y - rho * x) / sqrt(1 - rho^2);

[n, ~] = size(heatFlow);
dt = T / n;

t_vec = dt * (n:-1:1)';
x_vec = x_new - points_x;
y_vec = y_new - points_y;

G_mat = G(x_vec, y_vec, t_vec);

density = G(x_new, y_new, T) + dt / 2 * sum(sum(G_mat .* heatFlow .* S));

end