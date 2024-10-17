function DD = DiagonalDominante(A)
% DD = DiagonalDominante(A)
% Calcula si la matriz A es diagonal dominante, 1 si lo es, 0 si no
% PARAMETROS:
% A -> Matriz A del sistema debe ser simetrica y definida positiva   

    [n, m] = size(A);
    i = 1; DD = 1;

    while DD == 1 && i <= n
        sum = 0;
        for j = 1:m
            if j ~= i
                sum = sum + abs(A(i, j));
            end
        end

        if abs(A(i, i)) < sum
            DD = 0;
        end
    end
end