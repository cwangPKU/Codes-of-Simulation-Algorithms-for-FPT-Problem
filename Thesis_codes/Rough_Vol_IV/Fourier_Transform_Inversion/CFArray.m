function result = CFArray(H, rho, nu, lambda, theta, V0, u, t, r, d, n)
    [hs, Fs] = AdamsFractionalEquationArray(H, rho, nu, lambda, u, t, n);
    % disp(size(hs));
    % disp(size(Fs));
    hs(:, 1) = hs(:, 1) / 2;
    hs(:, n+1) = hs(:, n+1) / 2;
    G = sum(hs, 2) * t / n;

    Fs(:, 1) = Fs(:, 1) / 2;
    Fs(:, n+1) = Fs(: ,n+1) / 2;
    H = sum(Fs, 2) * t / n;
    
    result = exp(((r - d) .* t .* u .* 1i) + (lambda * theta * G + V0 * H).');
end