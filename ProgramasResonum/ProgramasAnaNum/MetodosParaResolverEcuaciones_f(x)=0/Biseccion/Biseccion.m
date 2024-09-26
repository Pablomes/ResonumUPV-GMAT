function [sol, incr, iter, ACOC] = Biseccion(f, a, b, tol, maxiter)

% [sol, incr, iter, ACOC] = Biseccion(f, a, b, tol, maxiter)
% Aproxima la solucion f(x) = 0 usando el metodo de biseccion
% PARAMETROS:
% f -> funcion anonima a resolver
% [a, b] -> intervalo en el que se encuentra la solucion
% tol -> tolerancia. Detiene ejecucion cuando se alcanza un incremento
% menor
% maxiter -> numero de iteraciones tras las que cesa la ejecucion

x0 = 0;
x1 = a;
iter = 0;
incr = [1];
ACOC = [];

while and(incr(end) > tol, iter < maxiter)
    x0 = (a + b) / 2;
    
    if f(a) * f(x0) < 0
        b = x0;
    elseif f(b) * f(x0) < 0
        a = x0;
    end

    incr = [incr, abs(x1 - x0)];
    iter = iter + 1;
    x1 = x0;
end

if iter>= maxiter
    sol = 'No ha convergido';
else
    sol = x0;
    incr = incr(:, 2:end);
    ACOC = fACOC(incr);
    Tasa(incr);
end