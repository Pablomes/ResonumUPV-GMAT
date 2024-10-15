function [sol, incr, iter, ACOC] = SOR1(A, b, x0, w, tol, maxiter)
% [sol, incr, iter, ACOC] = SOR1(A, b, x0, w, tol, maxiter)
% Aproxima la solucion Ax = b usando el metodo de sobre-relajacion SOR1
% PARAMETROS:
% A -> Matriz A del sistema debe ser cuadrada y con diagonal sin ceros
% b -> Vector independiente del sistema
% x0 -> vector solucion de inicio
% w -> parametro de relajacion. 1 < w < 2. w = 1 representa Gauss-Seidel
% tol -> tolerancia. Detiene ejecucion cuando se alcanza un incremento
% menor
% maxiter -> numero de iteraciones tras las que cesa la ejecucion
    
    b = b(:); x0 = x0(:); n = size(A, 1);
    incr = zeros(maxiter + 1, 1); incr(1) = 1;
    iter = 0;
    ACOC = [];


    while iter < maxiter && incr(iter + 1) > tol
        x1 = x0;
        for i = 1:n
            sigma = 0;
            for j = 1:n
                if j ~= i
                    sigma = sigma + A(i, j) * x1(j);
                end
            end
            x1(i) = (1 - w) * x1(i) + w * (b(i) - sigma) / A(i, i);
        end

        iter = iter + 1;
        incr(iter + 1) = max(x1 - x0);%norm(x1 - x0, inf);
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