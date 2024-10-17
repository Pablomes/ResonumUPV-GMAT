function [A, b] = Practica2_4_matriz(n)
% [A, b] = Practica1_2matriz(n)
% Genera la matriz y el vector del problema 4 de la practica 2
% PARAMETROS:
% n -> Tamaño del vector y de la matriz cuadrada

    A = zeros(n, n); b = ones(n, 1);

    for i = 1:n
        for j = i:n
            if i == j
                A(i, j) = 50;
            else
                A(i, j) = (j - i) / (5 * n);
                A(j, i) = A(i, j);
            end
        end
    end

    b(2:2:end) = 0;
end