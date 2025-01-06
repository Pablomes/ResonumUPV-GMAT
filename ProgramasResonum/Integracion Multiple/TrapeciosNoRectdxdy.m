function I = TrapeciosNoRectdxdy(F, ay, by, c, d, n, m)
% I = TrapeciosNoRectdxdy(F, ay, by, c, d, n, m)
% Aproxima la integral sobre el recinto [a(y), b(y)] x [c, d] haciendo
% int(F(x, y) dx dy) por el metodo de trapecios
% PARAMETROS:
% F -> Funcion a integrar. Debe ser anonima del tipo F = @(x, y) ...
% [a(y), b(y)] -> Dominio de integracion de la x. Deben ser del tipo
% a = @(y) ... y b = @(y) ...
% [c, d] -> Dominio de integracion de la y
% n -> numero de subintervalos para la x
% m -> numero de subintervalos para la y

    k = (d - c) / m;

    yk = (c:k:d)';

    sz = max(n + 1, m + 1);

    pesos = 2 * ones(sz, 1);
    pesos(1) = 1;
    pesos(end) = 1;

    Ix = zeros(m + 1, 1);
    for i = 1:m + 1
        if ay(yk(i)) == by(yk(i))
            continue;
        end

        h = (by(yk(i)) - ay(yk(i))) / n;
        xk = (ay(yk(i)):h:by(yk(i)))';

        z = F(xk, yk(i));
        Ix(i) = (h/2) * sum(pesos(1:m + 1) .* z);
    end

    I = (k/2) * sum(pesos(1:n + 1) .* Ix);
end