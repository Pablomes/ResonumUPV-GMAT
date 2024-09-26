function [O, E] = OrdenConvergenciaHeunSistemasSinSolucion(F, t0, T, Y0, N, tol)
% [O, E] = OrdenConvergenciaHeunSistemasSinSolucion(F, t0, T, Y0, N, tol)
% Obtiene el orden de convergencia del metodo de sistemas de Heun sin solucion
% conocida
% PARAMETROS:
% F -> sistema de EDOs de la forma f(t, y(t)). Debe ser una función anonima de Matlab
%       de formato f = @ (t, y) [... ; ...]
% t0 -> t conidicion inicial
% T -> t final. Suele ser el punto del que se quiere la aproximacion
% Y0 -> lista de y(t0) condiciones iniciales
% N -> Numero de subintervalos usados. Cuanto mas mejor sera la
%       aproximacion
% tol -> tolerancia de parada. Detendrá el programa cuando el orden de
%       convergencia no varie mas que tol

E = []; O = []; err = tol + 1; n = N; YK1 = []; YK2 = [];

[Yk1, ~] = MetodoHeunSistemas(F, t0, T, Y0, n);
n = n*2;
[Yk2, ~] = MetodoHeunSistemas(F, t0, T, Y0, n);
temp = Yk2(:, 1:2:end);
E = [E, sum(max(abs(Yk1 - temp), [], 2))];
Yk1 = Yk2;
n = n * 2;

while err > tol
    [Yk2, ~] = MetodoHeunSistemas(F, t0, T, Y0, n);
    temp = Yk2(:, 1:2:end);
    E = [E, sum(max(abs(Yk1 - temp), [], 2))];
    Yk1 = Yk2;

    O = [O, log2(E(:, end - 1) ./ E(:, end))];
    if length(O) > 1
        err = abs(O(end) - O(end - 1));
    end

    n = n*2;
end

end