function I = MonteCarloNoRectdydx(F, a, b, cx, dx, n)
% I = MonteCarloNoRectdydx(F, a, b, cx, dx, n)
% Aproxima la integral sobre el rectangulo [a, b] x [c(x), d(x)] haciendo
% int(F(x, y) dy dx) por el metodo de MonteCarlo
% PARAMETROS:
% F -> Funcion a integrar. Debe ser anonima del tipo F = @(x, y) ...
% [a, b] -> Dominio de integracion de la x
% [c(x), d(x)] -> Dominio de integracion de la y. Deben ser del tipo 
% c = @(x) ... y d = @(x) ...
% n -> numero de puntos aleatorios

    xk = a + (b - a) * rand(n, 1);
    yk = rand(n, 1);
    
    for i = 1:n
        yk(i) = yk(i) * (dx(xk(i)) - cx(xk(i))) + cx(xk(i));
    end

    sum = 0;

    for i = 1:n
        sum = sum + F(xk(i), yk(i));
    end

    I = (integral(dx, a, b) - integral(cx, a, b)) * (1 / n) * sum;
end