function [sol, incr, iter, ACOC] = Jacobi(A, b, x0, tol, maxiter)
% [sol, incr, iter, ACOC] = Jacobi(A, b, x0, tol, maxiter)
% Aproxima la solucion Ax = b usando el metodo estacionario de Jacobi
% PARAMETROS:
% A -> Matriz A del sistema. Debe ser cuadrada y con diagonal sin ceros
% b -> Vector independiente del sistema
% x0 -> vector solucion de inicio
% tol -> tolerancia. Detiene ejecucion cuando se alcanza un incremento
% menor
% maxiter -> numero de iteraciones tras las que cesa la ejecucion


    b = b(:); x0 = x0(:); n = size(A, 1);
    iter = 0;
    incr = zeros(maxiter + 1, 1); incr(1) = 1;
    ACOC = [];
    Di = diag(1 ./ diag(A));
    T = eye(n) - Di * A;

    while and(incr(iter + 1) > tol, iter < maxiter)
        x1 = T * x0 + Di * b;
        
        iter = iter + 1;
        incr(iter + 1) = norm(x1 - x0);
        x0 = x1;
    end

    if iter >= maxiter
        sol = 'No ha convergido';
    else
        sol = x0;
        incr = incr(2:iter + 1);
        ACOC = fACOC(incr);
    end
end