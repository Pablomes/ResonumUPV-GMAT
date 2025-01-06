function I = SimpsonRect(F, a, b, c, d, n, m)
% I = SimpsonRect(F, a, b, c, d, n, m)
% Aproxima la integral sobre el rectangulo [a, b] x [c, d] haciendo
% int(F(x, y) dy dx) por el metodo de Simpson
% PARAMETROS:
% F -> Funcion a integrar. Debe ser anonima del tipo F = @(x, y) ...
% [a, b] -> Dominio de integracion de la x
% [c, d] -> Dominio de integracion de la y
% n -> mitad del numero de subintervalos para la x
% m -> mitad del numero de subintervalos para la y

    h = (b - a) / (2 * n);
    k = (d - c) / (2 * m);

    xk = (a:h:b)';
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

    Iy = zeros(2 * n + 1, 1);
    for i = 1:2 * n + 1
        z = F(xk(i), yk);
        Iy(i) = (k/3) * sum(pesos(1:2 * m + 1) .* z);
    end

    I = (h/3) * sum(pesos(1:2 * n + 1) .* Iy);
end