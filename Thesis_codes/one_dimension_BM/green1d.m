function ans_green = green1d(x, t, a, b, N)
    n = 1:N;
    [n_grid, x_grid] = meshgrid(n, x);
    term = 2/(b-a) * exp(-(n_grid*pi/(b-a)).^2 * t/2) ...
        .* sin(-a/(b-a)*n_grid*pi) ...
        .* sin((x_grid-a)/(b-a).*n_grid*pi);
    ans_green = sum(term, 2);
end