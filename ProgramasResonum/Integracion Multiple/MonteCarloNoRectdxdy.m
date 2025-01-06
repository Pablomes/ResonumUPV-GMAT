function I = MonteCarloNoRectdxdy(F, ay, by, c, d, n)
% I = MonteCarloNoRectdydx(F, a, b, cx, dx, n)
% Aproxima la integral sobre el rectangulo [a, b] x [c(x), d(x)] haciendo
% int(F(x, y) dy dx) por el metodo de MonteCarlo
% PARAMETROS:
% F -> Funcion a integrar. Debe ser anonima del tipo F = @(x, y) ...
% [a(y), b(y)] -> Dominio de integracion de la x. Deben ser del tipo
% a = @(y) ... y b = @(y) ...
% [c, d] -> Dominio de integracion de la y
% n -> numero de puntos aleatorios

    xk = rand(n, 1);
    yk = c + (d - c) * rand(n, 1);
    
    for i = 1:n
        xk(i) = xk(i) * (by(yk(i)) - ay(yk(i))) + ay(yk(i));
    end

    sum = 0;

    for i = 1:n
        sum = sum + F(xk(i), yk(i));
    end

    I = (integral(by, c, d) - integral(ay, c, d)) * (1 / n) * sum;
end