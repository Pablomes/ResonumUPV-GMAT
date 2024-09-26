function [sol,incr,iter,ACOC] = ChebyshevHalley (beta,f,df,d2f,x0,tol,maxiter)
% [sol, incr, iter, ACOC] = Chebyshevhalley(beta, f, df, d2f, x0, tol, maxiter)
% Aproxima la solucion f(x) = 0 usando el metodo de punto fijo de
% Chebyshev-Halley
% PARAMETROS:
% beta ->   si beta = 0     ==> Metodo de Chebyshev
%           si beta = 1/2   ==> metodo de Halley
%           si beta = 1     ==> metodo de super-Halley
%           si beta -> inf  ==> metodo de Newton
% f -> funcion anonima a resolver
% df -> primera derivada de la funcion
% d2f -> segunda derivada de la funcion
% x0 -> punto de inicio
% tol -> tolerancia. Detiene ejecucion cuando se alcanza un incremento
% menor
% maxiter -> numero de iteraciones tras las que cesa la ejecucion

% 1. Inicializar variables
incr = [1];
iter = 0;
ACOC = [];

L = @(x) (f(x) * d2f(x)) / (df(x).^2);
% 2. Bucle iterativo
while and(incr(end) > tol, iter < maxiter)
    x1 = x0 - (1 + 0.5 * (L(x0) / (1 - beta * L(x0)))) * (f(x0) / df(x0));
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