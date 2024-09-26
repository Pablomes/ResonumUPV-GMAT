function [yk, xk] = DiferenciaFinitaLinealO2(p, q, r, a, b, alpha, beta, N)
% [yk, xk] = DiferenciaFinitaLinealO2(p, q, r, a, b, alpha, beta, N)
% Obtiene la solucion del PC lineal usando el metodo de Diferencias Finitas para EDOS
% de orden 2 de la forma y''(x) = p(x)y'(x) + q(x)y(x) + r(x) con
% condiciones de Dirichlet.
% PARAMETROS:
% p -> coeficiente del término y'(x). Debe ser una funcion anonima
% q -> coeficiente del termino y(x). Debe ser una función anonima
% r -> coeficiente del termino independiente. Debe ser una funcion anonima.
% [a, b] -> Dominio del PC
% alpha -> valor de y(a)
% beta -> valor de y(b)
% N -> numero de nodos intermedios

h = (b - a) / (N + 1); xk = a:h:b;
DSA = zeros(1, N-1); DPA = zeros(1, N); DIA = zeros(1, N-1); d = zeros(1, N);
d(1) = (1 + h/2 * p(xk(2))) * alpha; d(N) = (1 - h/2 * p(xk(N+1))) * beta;

for i = 1:N
    DPA(i) = 2 + h^2 * q(xk(i + 1));
    d(i) = d(i) - h^2 * r(xk(i + 1));

    if i < N
        DSA(i) = -1 + h/2 * p(xk(i + 1));
        DIA(i) = -1 - h/2 * p(xk(i + 2));
    end
end

yk = ResolucionDirectaCrout(DPA(:), DIA(:), DSA(:), d(:));
yk = [alpha, yk', beta];

hold off
plot(xk, yk, "o--")
axis padded
hold on

end