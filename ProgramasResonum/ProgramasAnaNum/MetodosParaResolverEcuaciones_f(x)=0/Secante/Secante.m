function [sol,incr,iter,ACOC] = Secante (f,x0,x1,tol,maxiter)
% [sol, incr, iter, ACOC] = Secante(f, x0, x1, tol, maxiter)
% Aproxima la solucion f(x) = 0 usando el metodo se la secante
% PARAMETROS:
% f -> funcion anonima a resolver
% x0 -> primer punto de inicio
% x1 -> segundo punto de inicio
% tol -> tolerancia. Detiene ejecucion cuando se alcanza un incremento
% menor
% maxiter -> numero de iteraciones tras las que cesa la ejecucion

% 1. Inicializar variables
incr = [1];
iter = 0;
ACOC = [];
% 2. Bucle iterativo
while and(incr(end) > tol, iter < maxiter)
    x2 = x1 - ((f(x1) * (x1 - x0)) / (f(x1) - f(x0)));
    % actualizar valores para condición de parada
    incr = [incr, abs(x2 - x1)];
    iter = iter + 1;
    x0 = x1; x1 = x2;
end
% 3. Procesar la solución
if iter>= maxiter
    sol = 'No ha convergido';
else
    sol = x1;
    incr = incr(:, 2:end);
    ACOC = fACOC(incr);
    Tasa(incr);
end