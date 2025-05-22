% This function serves to solve the fractional Riccati equation (4.7),
% which takes the general form of (4.13)

function [hs, Fs] = AdamsFractionalEquationArray(H, rho, nu, lambda, u, t, n)
% Here we allow u to be an array
    u = u.';

    hs = zeros(length(u), n+1);
    Fs = zeros(length(u), n+1);
    Fs(:, 1) = F(rho, nu, lambda, u, 0);
    delta = t / n;
    %disp(Fs);
    delta_Hpower = delta^(H+1/2);
    %gamma1 = gamma(H+1/2);
    gamma3 = gamma(H+3/2);
    gamma5 = gamma(H+5/2);

    %disp([delta_Hpower, gamma1, gamma3, gamma5]);
    

    for k = 0:n-1
        as = zeros(1, k+2);
        as(1) = delta_Hpower / gamma5 * (k^(H+3/2) - (k-H-1/2) * (k+1)^(H+1/2));
        %fa = @(j) delta_Hpower / gamma5 .* ((k-j+2).^(H+3/2) + (k-j).^(H+3/2) - 2*(k-j+1).^(H+3/2));
        %as(2:k+1) = fa(1:k);
        for j = 1:k
            as(j+1) = delta_Hpower / gamma5 * ((k-j+2)^(H+3/2) + (k-j)^(H+3/2) - 2*(k-j+1)^(H+3/2));
        end
        as(k+2) = delta_Hpower / gamma5;

        bs = zeros(1, k+1);
        for j = 0:k
            bs(j+1) = delta_Hpower / gamma3 * ((k-j+1)^(H+1/2) - (k-j)^(H+1/2));
        end

        %disp(as);
        %disp(bs);

        hs(:, k+2) = sum(Fs(:, 1:k+1) .* bs, 2);
        %disp(size(hs(:, k+2)));
        %disp(hs(:, k+2));
        %disp(size(F(rho, nu, lambda, u, hs(:, k+2))));
        %disp(F(rho, nu, lambda, u, hs(:, k+2)));
        %disp(size(hs(:, k+2)));
        Fs(:, k+2) = F(rho, nu, lambda, u, hs(:, k+2));
        %disp([hs(k+2), Fs(k+2)]);
        %disp(size(hs(:, k+2)));
        hs(:, k+2) = sum(Fs(:, 1:k+2) .* as, 2);
        %disp(size(hs(:, k+2)));
        Fs(:, k+2) = F(rho, nu, lambda, u, hs(:, k+2));
    end

end


function result = F(rho, nu, lambda, u, h)
    % Here we allow u to be an array
    % disp(size(u));
    % disp(size(- u .* (u + 1i) ./ 2));
    result = - u .* (u + 1i) ./ 2 + (u .* rho .* nu .* 1i - lambda) .* h + nu.^2 .* h.^2 ./ 2;

end