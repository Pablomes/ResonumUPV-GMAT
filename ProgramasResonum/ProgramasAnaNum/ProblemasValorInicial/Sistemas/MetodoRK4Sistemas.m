function [Yk, tk] = MetodoRK4Sistemas(F, t0, T, Y0, N)
% [Yk, tk] = MetodoRK4Sistemas(F, t0, T, Y0, N)
% Obtiene la solucion del PVI sistema usando el metodo de Runge-Kutta. Devuelve tambien
% los nodos usados.
% PARAMETROS:
% F -> sistema de EDOs de la forma f(t, y(t)). Debe ser una función anonima de Matlab
%       de formato f = @ (t, y) [... ; ...]
% t0 -> t conidicion inicial
% T -> t final. Suele ser el punto del que se quiere la aproximacion
% Y0 -> lista de y(t0) condiciones iniciales
% N -> Numero de subintervalos usados. Cuanto mas mejor sera la
%       aproximacion

h = (T - t0) / N; Yk = zeros(length(Y0), N + 1); tk = linspace(t0, T, N + 1);
Yk(:, 1) = Y0(:);

for i = 1:N
    K1 = F(tk(i), Yk(:, i));
    K2 = F(tk(i) + h/2, Yk(:, i) + h/2 * K1);
    K3 = F(tk(i) + h/2, Yk(:, i) + h/2 * K2);
    K4 = F(tk(i+1), Yk(:, i) + h * K3);

    Yk(:, i + 1) = Yk(:, i) + h / 6 * (K1 + 2 * K2 + 2 * K3 + K4);
end

for i = 1:length(Y0)
    plot(tk, Yk(i, :), "o--")
    hold on
end
axis padded

end