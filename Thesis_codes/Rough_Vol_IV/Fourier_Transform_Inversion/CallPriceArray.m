function C = CallPriceArray(S, k, H, rho, nu, lambda, theta, V0, t, r, d, n, left, right, N)
    [Pi1, Pi2] = PiArray(k, H, rho, nu, lambda, theta, V0, t, r, d, n, left, right, N);
    C = S .* Pi1 - exp(-r * t) * exp(k) .* Pi2;
end