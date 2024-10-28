function [Q, R] = TriangulizacionHouseholder(A)
% [Q, R] = TriangulizacionHouseholder(A)
% Utiliza las matrices de Householder para triangulizar la matriz
% PARAMETROS:
% A -> Matriz a triangulizar

    [m, n] = size(A);

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
end