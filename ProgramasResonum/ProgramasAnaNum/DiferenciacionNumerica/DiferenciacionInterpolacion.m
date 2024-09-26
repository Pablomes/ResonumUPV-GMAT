function [dp, dVal] = DiferenciacionInterpolacion(f, x0, n, h, tipo)
% [dp, dVal] = DiferenciacionInterpolacion(f, x0, n, h, tipo)
% Aproxima el diferencial a base de interpolacion usando la funcion
% PARAMETROS:
% f -> funcion a diferenciar, debe ser una expresion simbolica
% x0 -> punto en el que se quiere la derivada
% n -> grado del polinomio de aproximacion de la derivada
% h -> distancia entre los nodos
% tipo -> tipo de creacion de los nodos
%   1 -> Regresiva
%   2 -> Progresiva
%   Otro -> Central

syms x; nodos = 0; fi = []; D = [];

if tipo == 1
    nodos = x0 - (h * n) : h : x0;
elseif tipo == 2
    nodos = x0 : h : x0 + (h * n);
else
    nodos = x0 - (h * n) : 2 * h : x0 + (h * n);
end

nodos = nodos(:);

fi = subs(f, x, nodos);
D = sym([fi, zeros(n + 1, n)]); b = sym([fi(1);zeros(n+1, 1)]);
p = fi(1); mult = 1;

for col = 2:(n + 1)
    for row = col:(n + 1)
        D(row, col) = (D(row, col-1) - D(row-1, col-1))/(nodos(row) - nodos(row - col + 1));
    end
    b(col) = D(col, col);
    mult = mult * (x - nodos(col - 1));
    p = p + (b(col) * mult);
end

dp = diff(p);
dVal = vpa(subs(dp, x, x0));

end