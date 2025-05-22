function output = G(x1, x2, t)


output = exp(-1 * (x1.^2 + x2.^2) ./ (2*t)) ./ (2 * pi * t);

end