function [sol, incr, iter, ACOC] = GaussSeidel(A, b, x0, tol, maxiter)
% [sol, incr, iter, ACOC] = GaussSeidel(A, b, x0, tol, maxiter)
% Aproxima la solucion Ax = b usando el metodo estacionario de Gauss-Seidel
% PARAMETROS:
% A -> Matriz A del sistema debe ser cuadrada y con diagonal sin ceros
% b -> Vector independiente del sistema
% x0 -> vector solucion de inicio
% tol -> tolerancia. Detiene ejecucion cuando se alcanza un incremento
% menor
% maxiter -> numero de iteraciones tras las que cesa la ejecucion

    b = b(:); x0 = x0(:); n = size(A, 1);
    iter = 0;
    incr = zeros(maxiter + 1, 1); incr(1) = 1;
    ACOC = [];
    DL = tril(A); U = triu(A, 1);

    while and(incr(iter + 1) > tol, iter < maxiter)
        % Sustitucion directa
        x1 = zeros(n, 1);
        z = b - U * x0;
        x1(1) = z(1) / DL(1, 1);

        for i = 2:n
            x1(i) = (z(i) - DL(i, 1:i - 1) * x1(1:i - 1)) / DL(i, i);
        end
        
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