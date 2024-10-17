function [sol, incr, iter, ACOC] = GradienteConjugado(A, b, x0, tol, maxiter)
% [sol, incr, iter, ACOC] = GradienteConjugado(A, b, x0, tol, maxiter)
% Aproxima la solucion Ax = b usando el metodo del gradiente conjugado
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
    r0 = b - A * x0;
    d0 = r0; t0 = (r0' * r0) / (r0' * A * r0);
    x1 = x0 + t0 * d0; r1 = b - A * x1;
    incr(1) = sqrt(r1' * r1);
    d1 = r1 - ((d0' * A * r1) / (d0' * A * d0)) * d0;

    while iter < maxiter && incr(iter + 1) > tol
        t = (r1' * d1) / (d1' * A * d1);
        x = x1 + t * d1;
        r = b - A * x;
        d = r - ((d1' * A * r) / (d1' * A * d1)) * d1;

        iter = iter + 1;
        incr(iter + 1) = sqrt(r' * r);
        d1 = d; x1 = x; r1 = r;
    end

    if iter >= maxiter
        sol = 'No ha convergido';
    else
        sol = x1;
        incr = incr(1:iter + 1);
        ACOC = fACOC(incr);
    end

end