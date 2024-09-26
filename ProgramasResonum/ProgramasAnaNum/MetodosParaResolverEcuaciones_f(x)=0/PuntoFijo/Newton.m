function [sol,incr,iter,ACOC] = Newton (f,df,x0,tol,maxiter)
% [sol, incr, iter, ACOC] = Newton(f, df, x0, tol, maxiter)
% Aproxima la solucion f(x) = 0 usando el metodo de punto fijo de Newton
% PARAMETROS:
% f -> funcion anonima a resolver
% df -> primera derivada de la funcion
% x0 -> punto de inicio
% tol -> tolerancia. Detiene ejecucion cuando se alcanza un incremento
% menor
% maxiter -> numero de iteraciones tras las que cesa la ejecucion

% 1. Inicializar variables
iter = 0;
incr = [1];
ACOC =  [];
% 2. Bucle iterativo
while and(incr(end) > tol, iter < maxiter)
    x1 = x0 - (f(x0) / df(x0)); % cambiar método por su expresión
    % actualizar valores para condición de parada
    incr = [incr, abs(x1 - x0)];
    iter = iter + 1;
    x0 = x1;
end
% 3. Procesar la solución
if iter>= maxiter
    sol = 'No ha convergido';
else
    sol = x0;
    incr = incr(:, 2:end);
    ACOC = fACOC(incr);
    Tasa(incr);
end