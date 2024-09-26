function [A, b] = Practica1_2matriz(n)
% [A, b] = Practica1_2matriz(n)
% Genera la matriz y el vector del problema 2 de la practica 1
% PARAMETROS:
% n -> Tamaño del vector y de la matriz cuadrada

    A = zeros(n, n); b = ones(n, 1);

    A = diag(10 * (1:n));
    A = triu(A + (1:n)' - (1:n) + 1);
    A = A + A' - (diag(diag(A))) - eye(n);

    b(2:2:end) = -1;

end