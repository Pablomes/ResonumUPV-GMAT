function [O, E] = OrdenConvergenciaEulerSistemasSolucion(F, t0, T, Y0, N, Y, tol)
% [O, E] = OrdenConvergenciaEulerSistemasSolucion(F, t0, T, Y0, N, Y, tol)
% Obtiene el orden de convergencia del metodo de sistemas de Euler con solucion
% conocida
% PARAMETROS:
% F -> sistema de EDOs de la forma f(t, y(t)). Debe ser una función anonima de Matlab
%       de formato f = @ (t, y) [... ; ...]
% t0 -> t conidicion inicial
% T -> t final. Suele ser el punto del que se quiere la aproximacion
% Y0 -> lista de y(t0) condiciones iniciales
% N -> Numero de subintervalos usados. Cuanto mas mejor sera la
%       aproximacion
% Y -> funciones solucion. Debe ser una lista de funciones anonimas de Matlab de formato
%       y = @(t) [... ; ...]
% tol -> tolerancia de parada. Detendrá el programa cuando el orden de
%       convergencia no varie mas que tol

E = []; O = []; err = tol + 1; n = N;

[Yk, tk] = MetodoEulerSistemas(F, t0, T, Y0, n);
E = [E, sum(max(abs(Y(tk) - Yk), [], 2))];
n = n * 2;

while err > tol
    [Yk, tk] = MetodoEulerSistemas(F, t0, T, Y0, n);
    E = [E, sum(max(abs(Y(tk) - Yk), [], 2))];
    O = [O, log2(E(:, end - 1) ./ E(:, end))];
    if length(O) > 1
        err = abs(O(end) - O(end - 1));
    end

    n = n*2;
end

end