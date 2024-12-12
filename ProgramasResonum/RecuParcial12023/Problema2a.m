function [Q, R] = Problema2a(A)

    [m, n] = size(A);

    Q = eye(m);
    R = A;

    for i = 1:n

        if (i>=(m - 1)) break; end

        v = R(i + 1:m, i);
        v(1) = v(1) + sqrt(v' * v);
        v
        input("Enter...");
        h = eye(m - i) - (2 / (v' * v)) * (v * v')
        input("Enter...");
        H = [eye(i), zeros(i, m - i); zeros(m - i, i), h]
        input("Enter...");

        R = H * R * H;
        Q = Q * H';
    end
end