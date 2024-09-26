function [p, coef] = SplineLineal(xi, fi)
% [p, coef] = SplineLineal(xi, fi)
% Aproxima los puntos usando un spline lineal
% PARAMETROS:
% xi -> lista de coordenadas x
% fi -> lista de coordenadas y

xi = xi(:); fi = fi(:); n = length(xi);
p = sym(zeros(n-1, 1)); syms x;

border = ((xi(n) - xi(1)) / n) / 4;
scatter(xi, fi)
hold on
coef = [];

for i = 1:(n-1)
    coef = [coef; [(fi(i+1) - fi(i)) / (xi(i+1) - xi(i)), fi(i)]];
    p(i) = fi(i) + ((fi(i+1) - fi(i)) / (xi(i+1) - xi(i))) * (x - xi(i));
    fplot(p(i), [xi(i), xi(i+1)])
    hold on
end
end