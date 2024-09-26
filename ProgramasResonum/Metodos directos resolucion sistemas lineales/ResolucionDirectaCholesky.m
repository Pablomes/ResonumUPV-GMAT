function sol = ResolucionDirectaCholesky(A, b)
% sol = ResolucionDirectaCholesky(A, b)
% Utiliza la factorizacion de Cholesky para resolver el sistema Ax = b
% PARAMETROS:
% A -> Matriz del sistema. Debe ser simetrica y definida positiva
% b -> Vector independiente del sistema

    % ENCONTRAMOS L MEDIANTE CHOLESKY

    n = size(A, 2);
    b = b(:);

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

    % SUSTITUCION DIRECTA Lz = b

    z = zeros(n, 1);

    z(1) = b(1) / L(1, 1);

    for i = 2:n
        z(i) = (b(i) - L(i, 1:i - 1) * z(1:i - 1)) / L(i, i);
    end

    % SUSTITUCION INVERSA L'x = z

    sol = zeros(n, 1);

    sol(n) = z(n) / L(n, n);

    for i = (n - 1):-1:1
        sol(i) = (z(i) - L(i + 1:n, i)' * sol(i + 1:n)) / L(i, i);
    end

end