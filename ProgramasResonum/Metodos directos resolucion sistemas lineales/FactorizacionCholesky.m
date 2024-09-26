function L = FactorizacionCholesky(A)
% L = FactorizacionCholesky(A)
% Obtiene la matriz L de la factorizacion de Cholesky A = LL'
% PARAMETROS:
% A -> Matriz a factorizar. Debe ser simetrica y definida positiva

    n = size(A, 2);
    L = zeros(n, n);

    L(1, 1) = sqrt(A(1, 1));

    L(:, 1) = A(:, 1) ./ L(1, 1);

    for i = 2:(n)
        L(i, i) = sqrt(A(i, i) - L(i, 1:(i - 1)) * L(i, 1:(i - 1))');
        
        if i < n
            for j = (i + 1):n
                L(j, i) = (A(i, j) - L(i, 1:(i - 1)) * L(j, 1:(i - 1))') / L(i, i);
            end
        end
    end
end