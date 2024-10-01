function sol = ResolucionDirectaHouseholder(A, b)
% sol = ResolucionDirectaHouseholder(A, b)
% Utiliza las matrices de Householder para resolver un sistema compatible
% determinado
% PARAMETROS:
% A -> Matriz del sistema
% b-> vector independiente del sistema

    [m, n] = size(A);
    b = b(:);

    % Triangulizacion Householder
    Q = eye(m);
    R = A;

    for i = 1:n

        if i>=m break; end

        v = R(i:m, i);
        v(1) = v(1) + sqrt(v' * v);
        h = eye(m + 1 - i) - (2 / (v' * v)) * (v * v');
        H = [eye(i - 1), zeros(i - 1, m + 1 - i); zeros(m + 1 - i, i - 1), h];

        R = H * R;
        Q = Q * H';
    end

    % Encontrar z

    z = Q' * b;

    % Sustitucion inversa

    sol = zeros(n, 1);

    sol(n) = z(n) / R(n, n);

    for i = (n - 1):-1:1
        sol(i) = (z(i) - R(i, i + 1:n) * sol(i + 1:n)) / R(i, i);
    end

end