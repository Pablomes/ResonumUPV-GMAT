function [dp, dVal] = DiferenciacionInterpolacionNodos(xi, fi, x0)
% [dp, dVal] = DiferenciacionInterpolacion(xi, fi, x0)
% Aproxima el diferencial a base de interpolacion usando nodos
% PARAMETROS:
% xi -> Componentes x de los nodos en orden
% fi -> Componentes y de los nodos en orden
% x0 -> punto en el que se quiere la derivada

syms x; D = [];
xi = xi(:); fi = fi(:); n = length(xi) - 1;

D = sym([fi, zeros(n + 1, n)]); b = sym([fi(1);zeros(n+1, 1)]);
p = fi(1); mult = 1;

for col = 2:(n + 1)
    for row = col:(n + 1)
        D(row, col) = (D(row, col-1) - D(row-1, col-1))/(xi(row) - xi(row - col + 1));
    end
    b(col) = D(col, col);
    mult = mult * (x - xi(col - 1));
    p = p + (b(col) * mult);
end

dp = diff(p);
dVal = vpa(subs(dp, x, x0));

end