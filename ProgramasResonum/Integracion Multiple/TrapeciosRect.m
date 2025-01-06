function I = TrapeciosRect(F, a, b, c, d, n, m)
% I = TrapeciosRect(F, a, b, c, d, n, m)
% Aproxima la integral sobre el rectangulo [a, b] x [c, d] haciendo
% int(F(x, y) dy dx) por el metodo de trapecios
% PARAMETROS:
% F -> Funcion a integrar. Debe ser anonima del tipo F = @(x, y) ...
% [a, b] -> Dominio de integracion de la x
% [c, d] -> Dominio de integracion de la y
% n -> numero de subintervalos para la x
% m -> numero de subintervalos para la y

    h = (b - a) / n;
    k = (d - c) / m;

    xk = (a:h:b)';
    yk = (c:k:d)';

    sz = max(n + 1, m + 1);

    pesos = 2 * ones(sz, 1);
    pesos(1) = 1;
    pesos(end) = 1;

    Iy = zeros(n + 1, 1);
    for i = 1:n + 1
        z = F(xk(i), yk);
        Iy(i) = (k/2) * sum(pesos(1:m + 1) .* z);
    end

    I = (h/2) * sum(pesos(1:n + 1) .* Iy);
end