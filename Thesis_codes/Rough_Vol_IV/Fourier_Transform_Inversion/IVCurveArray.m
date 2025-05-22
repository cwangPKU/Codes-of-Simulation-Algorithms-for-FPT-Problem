function output = IVCurveArray(logStrikes, H, rho, nu, lambda, theta, V0, t, r, d, n, left, right, N)
    S = ones(size(logStrikes));
    callPrices = CallPriceArray(S, logStrikes, H, rho, nu, lambda, theta, V0, t, r, d, n, left, right, N);
    strikes = exp(logStrikes);
    IVfun = @(n) IV(1, strikes(n), t, r, callPrices(n));
    output = arrayfun(IVfun, (1:length(logStrikes)));
end


function output = IV(S, K, T, r, C)
% calculate IV for call price C, using binary search    
    sigma1 = 1e-8;
    sigma2 = 10;
    if callPrice(S, K, T, r, sigma1) > C
        output = sigma1;
        %disp([callPrice(S, K, T, r, sigma1), C]);
    elseif callPrice(S, K, T, r, sigma2) < C
        output = sigma2;
    %if or(callPrice(S, K, T, r, sigma1) > C, callPrice(S, K, T, r, sigma2) < C)
    %    output = [];
        %disp([callPrice(S, K, T, r, sigma1), callPrice(S, K, T, r, sigma2), C]);
    else
        while sigma2 - sigma1 > 1e-6
            sigma = (sigma1 + sigma2) / 2;
            if callPrice(S, K, T, r, sigma) > C
                sigma2 = sigma;
            else
                sigma1 = sigma;
            end
        end
        output = (sigma1 + sigma2) / 2;
    end
end


function output = callPrice(S, K, T, r, sigma)
    d1 = (log(S / K) + (r + sigma^2 / 2) * T) / (sigma * sqrt(T));
    d2 = (log(S / K) + (r - sigma^2 / 2) * T) / (sigma * sqrt(T));
    output = S * normcdf(d1) - exp(-r * T) * K * normcdf(d2);
end