function [sol, incr, iter, ACOC] = GradienteII(A, b, x0, tol, maxiter)
% [sol, incr, iter, ACOC] = GradienteII(A, b, x0, tol, maxiter)
% Aproxima la solucion Ax = b usando el metodo del gradiente
% PARAMETROS:
% A -> Matriz A del sistema debe ser simetrica y definida positiva
% b -> Vector independiente del sistema
% x0 -> vector solucion de inicio
% tol -> tolerancia. Detiene ejecucion cuando se alcanza un incremento
% menor
% maxiter -> numero de iteraciones tras las que cesa la ejecucion
    

    b = b(:); x0 = x0(:); n = size(A, 1);
    incr = zeros(maxiter + 1, 1);
    iter = 0;
    ACOC = [];
    r0 = b - A * x0; incr(1) = sqrt(r0' * r0);

    while iter < maxiter && incr(iter + 1) > tol
        v = A * r0;
        t = (r0' * r0) / (r0' * v);

        x1 = x0 + t * r0;
        r1 = r0 - t * v;

        iter = iter + 1;
        incr(iter + 1) = sqrt(r1' * r1);
        x0 = x1; r0 = r1;
    end

    if iter >= maxiter
        sol = 'No ha convergido';
    else
        sol = x0;
        ACOC = fACOC(incr);
    end
end