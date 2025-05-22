function [prob_a, prob_b, prob_in] = prob1d(t, a, b ,N)
    k = (-N:N+1)';
    prob_a = sum(normcdf(((2*k + 2) * (b - a) + a) ./ sqrt(t)) - normcdf((2 * k * (b - a) - a) ./ sqrt(t)), 1);
    prob_b = sum(normcdf(((2*k + 1) * (b - a) - a) ./ sqrt(t)) - normcdf(((2 * k + 1) * (b - a) + a) ./ sqrt(t)), 1);
    prob_in = 1 - prob_a - prob_b;
end