function sol = ResolucionDirectaFactCrout(A, b)
% sol = ResolucionDirectaFactCrout(A, b)
% Utiliza la factorizacion de Crout para resolver el sistema Ax = b con A
% simetrica
% PARAMETROS:
% A -> Matriz del sistema. Debe ser simetrica e invertible
% b -> Vector independiente del sistema

    n = size(A, 2);
    L = zeros(n, n);
    D = zeros(n, 1);

    % FACTORIZAMOS CROUT
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

    % SUSTITUCION DIRECTA Lz = b

    b = b(:);
    z = zeros(n, 1);

    z(1) = b(1);

    for i = 2:n
        z(i) = b(i) - L(i, 1:i - 1) * z(1:i - 1);
    end

    % DIAGONAL

    y = z ./ D;

    % SUSTITUCION INVERSA 

    y = y(:);
    sol = zeros(n, 1);
    U = L';

    sol(n) = y(n) / U(n, n);

    for i = (n - 1):-1:1
        sol(i) = (y(i) - U(i, i + 1:n) * sol(i + 1:n)) / U(i, i);
    end
end