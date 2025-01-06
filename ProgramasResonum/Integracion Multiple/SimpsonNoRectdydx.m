function I = SimpsonNoRectdydx(F, a, b, cx, dx, n, m)
% I = SimpsonNoRectdydx(F, a, b, cx, dx, n, m)
% Aproxima la integral sobre el recinto [a, b] x [c(x), d(x)] haciendo
% int(F(x, y) dy dx) por el metodo de Simpson
% PARAMETROS:
% F -> Funcion a integrar. Debe ser anonima del tipo F = @(x, y) ...
% [a, b] -> Dominio de integracion de la x
% [c(x), d(x)] -> Dominio de integracion de la y. Deben ser del tipo 
% c = @(x) ... y d = @(x) ...
% n -> mitad del numero de subintervalos para la x
% m -> mitad del numero de subintervalos para la y

    h = (b - a) / (2 * n);

    xk = (a:h:b)';

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
        if cx(xk(i)) == dx(xk(i))
            continue;
        end

        k = (dx(xk(i)) - cx(xk(i))) / (2 * m);
        yk = (cx(xk(i)):k:dx(xk(i)))';

        z = F(xk(i), yk);
        Iy(i) = (k/3) * sum(pesos(1:2 * m + 1) .* z);
    end

    I = (h/3) * sum(pesos(1:2 * n + 1) .* Iy);
end