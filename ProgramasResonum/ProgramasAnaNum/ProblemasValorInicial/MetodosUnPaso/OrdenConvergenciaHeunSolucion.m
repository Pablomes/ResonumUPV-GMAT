function [O, E] = OrdenConvergenciaHeunSolucion(f, t0, T, y0, N, y, tol)
% [O, E] = OrdenConvergenciaHeunSolucion(f, t0, T, y0, N, y, tol)
% Obtiene el orden de convergencia del metodo de Heun con solucion
% conocida
% PARAMETROS:
% f -> EDO de la forma f(t, y(t)). Debe ser una función anonima de Matlab
%       de formato f = @ (t, y) ...
% t0 -> t conidicion inicial
% T -> t final.
% y0 -> y(t0) condicion inicial 
% N -> Numero de subintervalos usados inicial.
% y -> funcion solucion. Debe ser una funcion anonima de Matlab de formato
%       y = @(t)
% tol -> tolerancia de parada. Detendrá el programa cuando el orden de
%       convergencia no varie mas que tol

E = []; O = []; err = tol + 1; n = N;

[yk, tk] = MetodoHeun(f, t0, T, y0, n);
E = [E, max(abs(y(tk) - yk))];
n = n * 2;

while err > tol
    [yk, tk] = MetodoHeun(f, t0, T, y0, n);
    E = [E, max(abs(y(tk) - yk))];

    O = [O, log2(E(end - 1) / E(end))];
    if length(O) > 1
        err = abs(O(end) - O(end - 1));
    end

    n = n*2;
end

end