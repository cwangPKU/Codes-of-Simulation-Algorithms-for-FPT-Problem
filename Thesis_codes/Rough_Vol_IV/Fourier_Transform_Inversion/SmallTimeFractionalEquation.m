function result = SmallTimeFractionalEquation(H, rho, nu, ~, u, t, n)

    gammas = gamma(1 + (H+1/2) .* (0:n));
    
    betas = zeros(1, n);
    betas(1) = -u * (u + 1i) / 2;
    for k = 1:n-1
        betas(k+1) = 1i * rho * u * gammas(k) / gammas(k+1) * betas(k);
        for a = 0:n-2
            b = n-2-a;
            new_item = betas(a+1) * betas(b+1) * gammas(a+1) * gammas(b+1) / (gammas(a+2) * gammas(b+2));
            betas(k+1) = betas(k+1) + new_item / 2;
        end
    end
    %disp(t.^(1:n).*(H+1/2));
    %disp(gammas);
    %disp(betas);
    %disp(gammas(1:n) ./ gammas(2:n+1) .* betas .* nu.^(1:n) .* t.^((1:n).*(H+1/2)));
    result = sum(gammas(1:n) ./ gammas(2:n+1) .* betas .* nu.^(0:n-1) .* t.^((1:n).*(H+1/2)));
    %disp(result);
end