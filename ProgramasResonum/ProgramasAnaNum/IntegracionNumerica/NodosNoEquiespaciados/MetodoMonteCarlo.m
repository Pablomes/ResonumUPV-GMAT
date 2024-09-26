function I = MetodoMonteCarlo(f, a, b, n)
% I = MetodoMonteCarlo(f, a, b, n)
% Obtiene la aproximacion de la integral mediante el metodo de Monte-Carlo
% en [a, b]
% PARAMETROS:
% f -> funcion a integrar. Debe ser una expresion simbolica.
% [a, b] -> Dominio de integracion
% n -> numero de puntos aleatorios generados. Cuanto mayor mejor sera la
%       aproximacion

syms x;
t = rand(1, n); xi = a + t * (b - a);

I = (b - a) / n * sum(subs(f, x, xi));
I = vpa(I, 7);

end