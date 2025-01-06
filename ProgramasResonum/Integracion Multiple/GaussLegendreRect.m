function I = GaussLegendreRect(F, a, b, c, d, n, m)
% I = GaussLegendreRect(F, a, b, c, d, n, m)
% Aproxima la integral sobre el rectangulo [a, b] x [c, d] haciendo
% int(F(x, y) dy dx) por la cuadratura de Gauss-Legendre
% PARAMETROS:
% F -> Funcion a integrar. Debe ser anonima del tipo F = @(x, y) ...
% [a, b] -> Dominio de integracion de la x
% [c, d] -> Dominio de integracion de la y
% n -> numero de nodos para x
% m -> numero de nodos para y

    [xNodos, xPesos] = NPLegendre(n);
    [yNodos, yPesos] = NPLegendre(m);

    xk = (a + b) / 2 + ((b - a) / 2) * xNodos;
    yk = (c + d) / 2 + ((d - c) / 2) * yNodos;
    
    Iy = zeros(n, 1);

    for i = 1:n
        z = F(xk(i), yk);
        Iy(i) = ((d - c) / 2) * sum(yPesos.*z);
    end

    I = ((b - a) / 2) * sum(xPesos.*Iy);
end