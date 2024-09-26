function [O, E] = OrdenConvergenciaEulerSinSolucion(f, t0, T, y0, N, tol)
% [O, E] = OrdenConvergenciaEulerSolucion(f, t0, T, y0, N, y, tol)
% Obtiene el orden de convergencia del metodo de Euler sin solucion
% conocida
% PARAMETROS:
% f -> EDO de la forma f(t, y(t)). Debe ser una función anonima de Matlab
%       de formato f = @ (t, y) ...
% t0 -> t conidicion inicial
% T -> t final.
% y0 -> y(t0) condicion inicial 
% N -> Numero de subintervalos usados inicial.
% tol -> tolerancia de parada. Detendrá el programa cuando el orden de
%       convergencia no varie mas que tol

E = []; O = []; err = tol + 1; n = N; yk1 = []; yK2 = [];

[yk1, ~] = MetodoEuler(f, t0, T, y0, n);
n = n * 2;
[yk2, ~] = MetodoEuler(f, t0, T, y0, n);
temp = yk2(1:2:end);
E = [E, max(abs(yk1 - temp))];
yk1 = yk2;
n = n * 2;

while err > tol
    [yk2, ~] = MetodoEuler(f, t0, T, y0, n);
    temp = yk2(1:2:end);
    E = [E, max(abs(yk1 - temp))];
    yk1 = yk2;

    O = [O, log2(E(end - 1) / E(end))];
    if length(O) > 1
        err = abs(O(end) - O(end - 1));
    end

    n = n*2;
end

end