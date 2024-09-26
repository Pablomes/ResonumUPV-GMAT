function [sol,incr,iter,ACOC] = PuntoFijo(g, x0, tol, maxiter)
% [sol, incr, iter, ACOC] = PuntoFijo(g, x0, tol, maxiter)
% Resuelve por el metodo de punto fijo
% PARAMETROS:
% g -> función tal que x = g(x)
% x0 -> x inicial
% tol -> tolerancia
% maxiter -> maximo numero de iteraciones


incr = [1];
iter = 0;
ACOC = [];

while and(incr(end) > tol, iter < maxiter)
    x1 = g(x0);
    % actualizar valores para condición de parada
    incr = [incr, abs(x1 - x0)];
    iter = iter + 1;
    x0 = x1;
end

if iter>= maxiter
    sol = 'No ha convergido';
else
    sol = x0;
    incr = incr(:, 2:end);
    ACOC = fACOC(incr);
    Tasa(incr);
end

end