function I = SimpsonCuboide(F, a, b, c, d, e, f, n, m, l)
% I = SimpsonRect(F, a, b, c, d, e, f, n, m, l)
% Aproxima la integral sobre el cuboide [a, b] x [c, d] x [e, f] haciendo
% int(F(x, y, z) dz dy dx) por el metodo de Simpson
% PARAMETROS:
% F -> Funcion a integrar. Debe ser anonima del tipo F = @(x, y, z) ...
% [a, b] -> Dominio de integracion de la x
% [c, d] -> Dominio de integracion de la y
% [e, f] -> Dominio de integracion de la z
% n -> mitad del numero de subintervalos para la x
% m -> mitad del numero de subintervalos para la y
% l -> mitad del numero de subintervalos para la y

    h = (b - a) / (2 * n);
    k = (d - c) / (2 * m);
    p = (f - e) / (2 * l);

    xk = (a:h:b)';
    yk = (c:k:d)';
    zk = (e:p:f)';

    sz = max(max(2 * n + 1, 2 * m + 1), 2 * l + 1);

    pesos = ones(sz, 1);
    for i = 2:sz - 1
        if mod(i, 2) == 0
            pesos(i) = 4;
        else
            pesos(i) = 2;
        end
    end
    
    Iy = zeros(2 * n + 1, 1);
    Iz = zeros(2 * m + 1, 1);
    for i = 1:2 * n + 1
        for j = 1: 2 * m + 1
            t = F(xk(i), yk(j), zk);
            Iz(j) = (p/3) * sum(pesos(1:2 * l + 1) .* t);
        end

        Iy(i) = (k/3) * sum(pesos(1:2 * m + 1) .* Iz);
    end

    I = (h/3) * sum(pesos(1:2 * n + 1) .* Iy);
end