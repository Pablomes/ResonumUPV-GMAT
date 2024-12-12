function [Q, R] = TridiagonalizacionHouseholder(A)
% [Q, R] = TridiagonalizacionHouseholder(A)
% Utiliza las matrices de Householder para tridiagonalizar la matriz
% PARAMETROS:
% A -> Matriz a tridiagonalizar. Debe ser simetrica

    [m, n] = size(A);

    Q = eye(m);
    R = A;

    for i = 1:n

        if (i>=(m - 1)) break; end

        v = R(i + 1:m, i);
        v(1) = v(1) + sqrt(v' * v);
        h = eye(m - i) - (2 / (v' * v)) * (v * v');
        H = [eye(i), zeros(i, m - i); zeros(m - i, i), h];

        R = H * R * H;
        Q = Q * H';
    end
end