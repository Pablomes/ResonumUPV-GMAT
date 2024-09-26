function [f, coef] = MinimosCuadradosNoLineales(xi, fi, tipo)
% [f, coef] = MinimosCuadradosNoLineales(xi, fi, tipo)
% Aproxima los puntos minimizando el error cuadratico en formas
% no lineales a base de linealizarlas primero
% PARAMETROS:
% xi -> lista de coordenadas x
% fi -> lista de coordenadas y
% tipo -> Numero de tipo
%TIPOS:
% 1: y = ce^(kx)
% 2: y = (a/x) + b
% 3: y = 1/(ax + b)^2 (NO FUNCIONA SIEMPRE)

xi = xi(:); fi = fi(:); syms x; 
n = length(xi);

border = ((xi(n) - xi(1)) / n) / 4;
scatter(xi, fi)
hold on

XI = zeros(length(xi), 1); FI = zeros(length(fi), 1);

if tipo == 1
    XI = xi;
    FI = log(fi);
elseif tipo == 2
    XI = 1./xi;
    FI = fi;
elseif tipo == 3
    XI = xi;
    FI = fi.^(-2);
else
    error("TIPO NO INCLUIDO.")
end

sumx = 0; sumxx = 0; sumy = 0; sumxy = 0;

for i = 1:n
    sumx = sumx + XI(i);
    sumxx = sumxx + (XI(i) * XI(i));
    sumy = sumy + FI(i);
    sumxy = sumxy + (XI(i) * FI(i));
end

coef = (1 / (n * sumxx - sumx * sumx)) * [sumxx * sumy - sumxy * sumx; n * sumxy - sumx * sumy];

if tipo == 1
    f = exp(coef(1)) * exp(coef(2) * x);
elseif tipo == 2
    f = coef(1) + coef(2) * (1/x);
elseif tipo == 3
    f = (coef(1) + coef(2) * x) ^ (-2);
end

fplot(f, [xi(1) - border, xi(n) + border])
hold on

end