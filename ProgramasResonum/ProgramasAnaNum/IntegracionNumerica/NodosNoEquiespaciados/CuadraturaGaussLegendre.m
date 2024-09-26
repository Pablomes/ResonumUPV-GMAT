function [xk, ck] = CuadraturaGaussLegendre(n)
% [xk, ck] = CuadraturaGaussLegendre(n)
% Obtiene los coeficientes c y las raices x de la cuadratura Gauss-Legendre
% en [-1, 1]
% PARAMETROS:
% n -> grado

syms x;
xk = zeros(n, 1); ck = zeros(n, 1);
p = LegendreP(n, x);
dp = diff(p);

xk = root(p, x);
for i = 1:n 
    %xk(i) = (1 - (1 / (8 * n^2)) + (1/(8 * n^3))) * cos((4 * i - 1) / (4 * n + 2) * pi);
    ck(i) = 2 / ((1 - xk(i)^2) * (subs(dp, x, xk(i)))^2);
end

xk = vpa(xk(end:-1:1), 5);
ck = ck(end:-1:1);

end