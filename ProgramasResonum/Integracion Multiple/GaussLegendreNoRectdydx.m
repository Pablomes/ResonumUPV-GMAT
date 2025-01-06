function I = GaussLegendreNoRectdydx(F, a, b, cx, dx, n, m)
% I = GaussLegendreNoRectdydx(F, a, b, cx, dx, n, m)
% Aproxima la integral sobre el recinto [a, b] x [c(x), d(x)] haciendo
% int(F(x, y) dy dx) por la cuadratura de Gauss-Legendre
% PARAMETROS:
% F -> Funcion a integrar. Debe ser anonima del tipo F = @(x, y) ...
% [a, b] -> Dominio de integracion de la x
% [c(x), d(x)] -> Dominio de integracion de la y. Deben ser del tipo 
% c = @(x) ... y d = @(x) ...
% n -> numero de nodos para x
% m -> numero de nodos para y

    [xNodos, xPesos] = NPLegendre(n);
    [yNodos, yPesos] = NPLegendre(m);

    xk = (a + b) / 2 + ((b - a) / 2) * xNodos;
    
    Iy = zeros(n, 1);

    for i = 1:n
        yk = (cx(xk(i)) + dx(xk(i))) / 2 + ((dx(xk(i)) - cx(xk(i))) / 2) * yNodos;

        z = F(xk(i), yk);
        Iy(i) = ((dx(xk(i)) - cx(xk(i))) / 2) * sum(yPesos.*z);
    end

    I = ((b - a) / 2) * sum(xPesos.*Iy);
end