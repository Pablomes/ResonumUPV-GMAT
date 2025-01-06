function I = TrapeciosNoRectdydx(F, a, b, cx, dx, n, m)
% I = TrapeciosNoRectdydx(F, a, b, cx, dx, n, m)
% Aproxima la integral sobre el recinto [a, b] x [c(x), d(x)] haciendo
% int(F(x, y) dy dx) por el metodo de trapecios
% PARAMETROS:
% F -> Funcion a integrar. Debe ser anonima del tipo F = @(x, y) ...
% [a, b] -> Dominio de integracion de la x
% [c(x), d(x)] -> Dominio de integracion de la y. Deben ser del tipo 
% c = @(x) ... y d = @(x) ...
% n -> numero de subintervalos para la x
% m -> numero de subintervalos para la y

    h = (b - a) / n;

    xk = (a:h:b)';

    sz = max(n + 1, m + 1);

    pesos = 2 * ones(sz, 1);
    pesos(1) = 1;
    pesos(end) = 1;

    Iy = zeros(n + 1, 1);
    for i = 1:n + 1
        if cx(xk(i)) == dx(xk(i))
            continue;
        end

        k = (dx(xk(i)) - cx(xk(i))) / m;
        yk = (cx(xk(i)):k:dx(xk(i)))';

        z = F(xk(i), yk);
        Iy(i) = (k/2) * sum(pesos(1:m + 1) .* z);
    end

    I = (h/2) * sum(pesos(1:n + 1) .* Iy);
end