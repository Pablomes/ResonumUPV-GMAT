function [I, xk, ck] = CuadraturaGaussLegendreAplicada(f, a, b, n)
% [I, xk, ck] = CuadraturaGaussLegendreAplicada(f, a, b, n)
% Obtiene la integral, los coeficientes c y las raices x de la cuadratura Gauss-Legendre
% en [a, b]
% PARAMETROS:
% f -> funcion a integrar. Debe ser una expresion simbolica
% [a, b] -> Dominio de integracion
% n -> grado

syms x t dx;
xk = zeros(n, 1); ck = zeros(n, 1);
p = LegendreP(n, t);
dp = diff(p);
f = f * dx;

tExp = ((b - a) * t + (b + a)) / 2;
dtExp = diff(tExp);
ft = subs(f, x, tExp);
ft = subs(ft, dx, dtExp);
I = 0;

xk = root(p, t);
for i = 1:n 
    %xk(i) = (1 - (1 / (8 * n^2)) + (1/(8 * n^3))) * cos((4 * i - 1) / (4 * n + 2) * pi);
    ck(i) = 2 / ((1 - xk(i)^2) * (subs(dp, t, xk(i)))^2);
    
    I = I + ck(i) * subs(ft, t, xk(i));
end

xk = vpa(xk(end:-1:1), 5);
ck = ck(end:-1:1);

end