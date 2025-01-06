function I = MonteCarloRect(F, a, b, c, d, n)
% I = MonteCarloRect(F, a, b, c, d, n)
% Aproxima la integral sobre el rectangulo [a, b] x [c, d] haciendo
% int(F(x, y) dy dx) por el metodo de MonteCarlo
% PARAMETROS:
% F -> Funcion a integrar. Debe ser anonima del tipo F = @(x, y) ...
% [a, b] -> Dominio de integracion de la x
% [c, d] -> Dominio de integracion de la y
% n -> numero de puntos aleatorios

    xk = a + (b - a) * rand(n, 1);
    yk = c + (d - c) * rand(n, 1);

    sum = 0;

    for i = 1:n
        sum = sum + F(xk(i), yk(i));
    end

    I = (b - a) * (d - c) * (1 / n) * sum;
end