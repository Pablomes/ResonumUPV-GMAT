function [f, coef] = MinimosCuadradosPolinomicos(xi, fi, grado)
% [f, coef] = MinimosCuadradosPolinomicos(xi, fi, grado)
% Aproxima los puntos minimizando el error cuadratico en un polinomio de
% grado dado
% PARAMETROS:
% xi -> lista de coordenadas x
% fi -> lista de coordenadas y
% grado -> grado del polinomio

xi = xi(:); fi = fi(:); syms x; 
n = length(xi);

border = ((xi(n) - xi(1)) / n) / 4;
scatter(xi, fi)
hold on

sumsx = zeros(1, 2*grado + 1);
sumsy = zeros(1, grado + 1);

sumsx(1) = n;

for i = 1:n
    for j = 0:2 * grado
        if j <= grado
            sumsy(j + 1) = sumsy(j + 1) + (fi(i) * (xi(i) ^ j));
        end
        
        if j > 0
            sumsx(j + 1) = sumsx(j + 1) + (xi(i) ^ j);
        end
    end
end

A = [];
b = [];

for i = 1:(grado + 1)
    A = [A; sumsx(i: i + grado)];
    b = [b; sumsy(i)];
end

coef = Cramer(A, b);

f = 0;

for i = 1:length(coef)
    f = f + (coef(i) * (x^(i-1)));
end

fplot(f, [xi(1) - border, xi(n) + border])
hold on

end