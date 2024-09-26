function [L, D] = FactorizacionCrout(A)
% L = FactorizacionCrout(A)
% Obtiene la matriz L, y el vector D de la factorizacion LDL^T donde D es la diagonal de la matriz
% diagonal de los pivotes
% PARAMETROS:
% A -> Matriz a factorizar. Debe ser simetrica e invertible, y permitir A=LU
    
    n = size(A, 2);
    L = zeros(n, n);
    D = zeros(n, 1);

    L(1, 1) = 1; D(1) = A(1, 1);

    for i = 2:n
        L(i, 1) = A(i, 1) / D(1);
        j = 2;
        while j < i
            L(i, j) = (A(i, j) - L(i, 1:i - 1) * diag(D(1:i - 1) * L(j, 1:i - 1))) / D(j);
            j = j + 1;
        end

        L(i, i) = 1;
        D(i) = A(i, i) - L(i, 1:i - 1) * diag(D(1:i - 1) * L(i, 1:i - 1));
    end
end