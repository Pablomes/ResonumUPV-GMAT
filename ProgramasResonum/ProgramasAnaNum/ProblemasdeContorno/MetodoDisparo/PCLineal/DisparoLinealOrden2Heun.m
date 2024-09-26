function [yk, xk] = DisparoLinealOrden2Heun(p, q, r, a, b, alpha, beta, N)
% [yk, xk] = DisparoLinealOrden2Heun(p, q, r, a, b, alpha, beta, N)
% Obtiene la solucion del PC usando el metodo de disparo lineal para EDOs
% de orden 2 de forma y''(x) = p(x)y'(x) + q(x)y(x) + r(x) con condiciones
% de Dirichlet usando el metodo de Heun
% PARAMETROS:
% p -> coeficiente del término y'(x). Debe ser una funcion anonima
% q -> coeficiente del termino y(x). Debe ser una función anonima
% r -> coeficiente del termino independiente. Debe ser una funcion anonima.
% [a, b] -> Dominio del PC
% alpha -> valor de y(a)
% beta -> valor de y(b)
% N -> numero de subintervalos

h = (b-a) / N; xk = a:h:b;

FU = @(x, U) [U(2); p(x) * U(2) + q(x) * U(1) + r(x)];
U0 = [alpha; 0];
[U, ~] = MetodoHeunSistemas(FU, a, b, U0, N);

FV = @(x, V) [V(2); p(x) * V(2) + q(x) * V(1)];
V0 = [0;1];
[V, ~] = MetodoHeunSistemas(FV, a, b, V0, N);

yk = U(1, :) + (beta - U(1, end))/V(1, end) * V(1, :);

hold off
plot(xk, yk, "o--")
axis padded
hold on

end