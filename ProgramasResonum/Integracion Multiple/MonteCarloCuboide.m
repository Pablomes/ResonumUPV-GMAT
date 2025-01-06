function I = MonteCarloCuboide(F, a, b, c, d, e, f, n)
% I = MonteCarloCuboide(F, a, b, c, d, e, f, n)
% Aproxima la integral sobre el rectangulo [a, b] x [c, d] x [e, f] haciendo
% int(F(x, y, z) dz dy dx) por el metodo de MonteCarlo
% PARAMETROS:
% F -> Funcion a integrar. Debe ser anonima del tipo F = @(x, y) ...
% [a, b] -> Dominio de integracion de la x
% [c, d] -> Dominio de integracion de la y
% [e, f] -> Dominio de integracion de la z
% n -> numero de puntos aleatorios

    xk = a + (b - a) * rand(n, 1);
    yk = c + (d - c) * rand(n, 1);
    zk = e + (f - e) * rand(n, 1);

    sum = 0;

    for i = 1:n
        sum = sum + F(xk(i), yk(i), zk(i));
    end

    I = (b - a) * (d - c) * (f - e) * (1 / n) * sum;
end