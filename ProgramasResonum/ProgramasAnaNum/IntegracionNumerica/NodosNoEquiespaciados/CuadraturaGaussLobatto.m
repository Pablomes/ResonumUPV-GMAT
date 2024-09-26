function [xk, ck] = CuadraturaGaussLobatto(n)
% [xk, ck] = CuadraturaGaussLobatto(n)
% Obtiene los coeficientes c y las raices x de la cuadratura Gauss-Lobatto
% en [-1, 1]
% PARAMETROS:
% n -> grado

syms x;
xk = zeros(n, 1); ck = zeros(n, 1);
p = LegendreP(n - 1, x);
dp = diff(p);

xk = root(dp, x);
for i = 1:n 
    
    if or (i == 1, i == n)
        ck(i) = 2/(n*(n-1));
    else
        ck(i) = 2 / ((n * (n - 1)) * (subs(p, x, xk(i - 1)) ^2));
    end
end

xk = vpa([-1; xk(end:-1:1); 1] , 5);
ck = ck(end:-1:1);

end