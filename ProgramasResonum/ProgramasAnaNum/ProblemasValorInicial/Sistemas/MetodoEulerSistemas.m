function [Yk, tk] = MetodoEulerSistemas(F, t0, T, Y0, N)
% [Yk, tk] = MetodoEulerSistemas(F, t0, T, Y0, N)
% Obtiene la solucion del PVI sistema usando el metodo de Euler. Devuelve tambien
% los nodos usados.
% PARAMETROS:
% F -> sistema de EDOs de la forma f(t, y(t)). Debe ser una función anonima de Matlab
%       de formato f = @ (t, y) [... ; ...]
% t0 -> t conidicion inicial
% T -> t final. Suele ser el punto del que se quiere la aproximacion
% Y0 -> lista de y(t0) condiciones iniciales
% N -> Numero de subintervalos usados. Cuanto mas mejor sera la
%       aproximacion

h = (T - t0) / N; Yk = zeros(length(Y0), N + 1); tk = zeros(1, N + 1);
Yk(:, 1) = Y0(:); tk(1) = t0;

for i = 1:N
    Yk(:, i + 1) = Yk(:, i) + h * F(tk(i), Yk(:, i));
    tk(i + 1) = tk(i) + h;
end

for i = 1:length(Y0)
    plot(tk, Yk(i, :), "o--")
    hold on
end
axis padded

end