function result = LargeTimeFractionalEquation(H, rho, nu, ~, u, t, n)

    Gammas = gamma(1 - (H+1/2) .* (0:n-1));
    %disp(Gammas);
    A = sqrt(u * (u+1i) - (rho^2) * (u^2));
    %disp(A);
    r = -1i * rho * u - A;

    gammas = zeros(1, n);
    gammas(1) = 1;
    gammas(2) = -1;
    %disp(gammas);
    for k = 2:n-1
        gammas(k+1) = -gammas(k) + r * Gammas(k+1) / (2 * A) * sum((gammas(2:k)./Gammas(2:k)) .* (gammas(k:-1:2)./Gammas(k:-1:2)));
    end
    %disp(gammas);
    %disp('hello');
    %disp(gammas .* nu.^(-1.*(0:n-1)) .* t.^(-1.*(0:n-1).*(H+1/2)));
    result = r * sum(gammas .* nu.^(-1.*(1:n)) .* t.^(-1.*(0:n-1).*(H+1/2)) ./ (A.^(0:n-1) .* Gammas));
end