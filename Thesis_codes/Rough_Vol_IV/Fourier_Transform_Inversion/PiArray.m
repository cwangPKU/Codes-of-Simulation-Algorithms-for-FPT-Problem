function [Pi1, Pi2] = PiArray(k, H, rho, nu, lambda, theta, V0, t, r, d, n, left, right, N)
    k = k.';
    % Here we allow k = log(K/S) be an array of shape (X, 1)
    f_minus1i = CFArray(H, rho, nu, lambda, theta, V0, -1i, t, r, d, n);
    u = linspace(left, right, N+1);

    cf2 = CFArray(H, rho, nu, lambda, theta, V0, u, t, r, d, n);
    integrand2 = real(exp(-1i * k * u ) .* (cf2 ./ (1i .* u)));
    integrand2(:, 1) = integrand2(:, 1) / 2;
    integrand2(:, N+1) = integrand2(:, N+1) / 2;
    integrand2 = integrand2 * (right-left) / N;
    Pi2 = sum(integrand2, 2) / pi + 1/2;
    Pi2 = Pi2.';
    
    cf1 = CFArray(H, rho, nu, lambda, theta, V0, u - 1i, t, r, d, n);
    integrand1 = real(exp(-1i * k * u) .* (cf1 ./ (1i .* u .* f_minus1i)));
    integrand1(:, 1) = integrand1(:, 1) / 2;
    integrand1(:, N+1) = integrand1(:, N+1) / 2;
    integrand1 = integrand1 * (right-left) / N;
    Pi1 = sum(integrand1, 2) / pi + 1/2;
    
    Pi1 = Pi1.';

    %disp(Pi1);
    %disp(Pi2);
    
    %hold off;
    %plot(u, integrand1);
    %hold on;
    %plot(u, integrand2);
    %disp([Pi1, Pi2]);
end