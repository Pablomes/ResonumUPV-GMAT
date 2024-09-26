function [I, xk, ck] = CuadraturaGaussLobattoAplicada(f, a, b, n)
% [I, xk, ck] = CuadraturaGaussLobattoAplicada(f, a, b, n)
% Obtiene la integral, los coeficientes c y las raices x de la cuadratura
% Gauss-Lobatto
% en [a, b]
% PARAMETROS:
% f -> funcion a integrar. Debe ser una expresion simbolica
% [a, b] -> Dominio de integracion
% n -> grado

syms x t dx;
xk = zeros(n, 1); ck = zeros(n, 1);
p = LegendreP(n - 1, x);
dp = diff(p);
f = f * dx;

tExp = ((b - a) * t + (b + a)) / 2;
dtExp = diff(tExp);
ft = subs(f, x, tExp);
ft = subs(ft, dx, dtExp);
I = 0;

xk = root(dp, x);
for i = 1:n 
    
    if or (i == 1, i == n)
        ck(i) = 2/(n*(n-1));
        if i == 1
            I = I + ck(i) * subs(ft, t, -1);
        else
            I = I + ck(i) * subs(ft, t, 1);
        end
    else
        ck(i) = 2 / ((n * (n - 1)) * (subs(p, x, xk(i - 1)) ^2));
        I = I + ck(i) * subs(ft, t, xk(i - 1));
    end
end

xk = vpa(xk(end:-1:1), 5);
ck = ck(end:-1:1);
I = vpa(I, 7);

end