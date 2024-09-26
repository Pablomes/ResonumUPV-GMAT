function [f, coef] = MinimosCuadradosLineales(xi, fi)
% [f, coef] = MinimosCuadradosLineales(xi, fi)
% Aproxima los puntos minimizando el error cuadratico en una
% recta
% PARAMETROS:
% xi -> lista de coordenadas x
% fi -> lista de coordenadas y

xi = xi(:); fi = fi(:); syms x; 
n = length(xi);

border = ((xi(n) - xi(1)) / n) / 4;
scatter(xi, fi)
hold on

sumx = 0; sumxx = 0; sumy = 0; sumxy = 0;

for i = 1:n
    sumx = sumx + xi(i);
    sumxx = sumxx + (xi(i) * xi(i));
    sumy = sumy + fi(i);
    sumxy = sumxy + (xi(i) * fi(i));
end

coef = (1 / (n * sumxx - sumx * sumx)) * [sumxx * sumy - sumxy * sumx; n * sumxy - sumx * sumy];

f = coef(1) + coef(2) * x;

fplot(f, [xi(1) - border, xi(n) + border])
hold on

end