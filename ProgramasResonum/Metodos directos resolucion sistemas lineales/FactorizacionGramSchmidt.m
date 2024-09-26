function [Q, R] = FactorizacionGramSchmidt(A)
% [Q, R] = FactorizacionGramSchmidt(A)
% Utiliza la factorizacion de vectores ortogonales de Gram-Schmidt para
% obtener las matrices Q y R
% PARAMETROS:
% A -> Matriz a factorizar

    [rows, cols] = size(A);

    V = zeros(rows, cols);
    Q = zeros(rows, cols);
    R = zeros(cols, cols);
    norms = zeros(cols, 1);
    prods = zeros(cols, 1);
    
    for i = 1:cols
        j = 1;
        ort = zeros(rows, 1);
        while j < i
            prods(j) = (A(:, i)' * V(:, j)) / norms(j);
            ort = ort + prods(j) * V(:, j);
            j = j + 1;
        end
        prods(i) = 1;

        V(:, i) = A(:, i) - ort;
        norms(i) = (V(:, i)' * V(:, i));
        Q(:, i) = V(:, i) / sqrt(norms(i));
        R(1:i, i) = prods(1:i) .* sqrt(norms(1:i));
    end
end