function ans_reflection = reflection1d(x, t, a, b, M, N)
    k = M:N;
    [k_grid, x_grid] = meshgrid(k, x);
    term1 = exp(-(x_grid - 2*k_grid*(b-a)).^2/(2*t));
    term2 = exp(-(x_grid - 2*b - 2*k_grid*(b-a)).^2/(2*t));
    ans_reflection = (term1 - term2)/sqrt(2*pi*t);
    ans_reflection = sum(ans_reflection, 2);
end