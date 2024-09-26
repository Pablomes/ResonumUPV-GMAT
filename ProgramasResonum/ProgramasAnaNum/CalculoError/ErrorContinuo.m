function [E2, Emax] = ErrorContinuo(f, p, a, b, h)
% [E2, Emax] = ErrorContinuo(f, p, a, b, h)
% Calcula el error cuadratico y el error maximo de una aproximacion de una
% funcion continua
% PARAMETROS:
% f -> funcion original
% p -> aproximacion de f
% [a, b] -> dominio de aproximacion
% h -> numero de trozos del dominio para Emax

    syms x;

    E2 = int((f - p) * (f - p), x, a, b);

    nodes = linspace(a, b, h + 1);
    n = length(nodes);

    Emax = 0;

    for i = 1:n
        e = abs(subs(f, x, nodes(i)) - subs(p, x, nodes(i)));

        if e > Emax
            Emax = e;
        end
    end
end