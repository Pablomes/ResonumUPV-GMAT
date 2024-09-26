function [yk, tk] = MetodoHeun(f, t0, T, y0, N)
% [yk, tk] = MetodoHeun(f, t0, T, y0, N)
% Obtiene la solucion del PVI usando el metodo de Heun. Devuelve tambien
% los nodos usados.
% PARAMETROS:
% f -> EDO de la forma f(t, y(t)). Debe ser una función anonima de Matlab
%       de formato f = @ (t, y) ...
% t0 -> t conidicion inicial
% T -> t final. Suele ser el punto del que se quiere la aproximacion
% y0 -> y(t0) condicion inicial 
% N -> Numero de subintervalos usados. Cuanto mas mejor sera la
%       aproximacion

h = (T - t0) / N; yk = zeros(N + 1, 1); tk = linspace(t0, T, N + 1); tk = tk(:);

yk(1) = y0;

for i = 1:N
    k1 = f(tk(i), yk(i));
    k2 = f(tk(i + 1), yk(i) + h * k1);

    yk(i + 1) = yk(i) + h/2 * (k1 + k2);
end

plot(tk, yk, "o--")
axis padded
hold on

end