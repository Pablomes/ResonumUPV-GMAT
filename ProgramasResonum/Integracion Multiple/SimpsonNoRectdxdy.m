function I = SimpsonNoRectdxdy(F, ay, by, c, d, n, m)
% I = SimpsonRect(F, a, b, c, d, n, m)
% Aproxima la integral sobre el recinto [a(y), b(y)] x [c, d] haciendo
% int(F(x, y) dy dx) por el metodo de Simpson
% PARAMETROS:
% F -> Funcion a integrar. Debe ser anonima del tipo F = @(x, y) ...
% [a(y), b(y)] -> Dominio de integracion de la x. Deben ser del tipo
% a = @(y) ... y b = @(y) ...
% [c, d] -> Dominio de integracion de la y
% n -> mitad del numero de subintervalos para la x
% m -> mitad del numero de subintervalos para la y

    k = (d - c) / (2 * m);

    yk = (c:k:d)';

    sz = max(2 * n + 1, 2 * m + 1);

    pesos = ones(sz, 1);
    for i = 2:sz - 1
        if mod(i, 2) == 0
            pesos(i) = 4;
        else
            pesos(i) = 2;
        end
    end

    Ix = zeros(2 * m + 1, 1);
    for i = 1:2 * m + 1
        if ay(yk(i)) == by(yk(i))
            continue;
        end

        h = (by(yk(i)) - ay(yk(i))) / (2 * n);
        xk = (ay(yk(i)):h:by(yk(i)));

        z = F(xk, yk(i));
        Ix(i) = (h/3) * sum(pesos(1:2 * m + 1) .* z);
    end

    I = (k/3) * sum(pesos(1:2 * n + 1) .* Ix);
end