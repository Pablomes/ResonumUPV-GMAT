function I = GaussLegendreNoRectdxdy(F, ay, by, c, d, n, m)
% I = GaussLegendreNoRectdxdy(F, ay, by, c, d, n, m)
% Aproxima la integral sobre el recinto [a(y), b(y)] x [c, d] haciendo
% int(F(x, y) dx dy) por la cuadratura de Gauss-Legendre
% PARAMETROS:
% F -> Funcion a integrar. Debe ser anonima del tipo F = @(x, y) ...
% [a(y), b(y)] -> Dominio de integracion de la x. Deben ser del tipo
% a = @(y) ... y b = @(y) ...
% [c, d] -> Dominio de integracion de la y
% n -> numero de nodos para x
% m -> numero de nodos para y

    [xNodos, xPesos] = NPLegendre(n);
    [yNodos, yPesos] = NPLegendre(m);

    yk = (c + d) / 2 + ((d - c) / 2) * yNodos;
    
    Ix = zeros(m, 1);

    for i = 1:m
        xk = (ay(yk(i)) + by(yk(i))) / 2 + ((by(yk(i))) - ay(yk(i)) / 2) * xNodos;

        z = F(xk, yk(i));
        Ix(i) = ((by(yk(i)) - ay(yk(i))) / 2) * sum(xPesos.*z);
    end

    I = ((d - c) / 2) * sum(yPesos.*Ix);
end