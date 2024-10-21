function [l, z, incr, iter, ACOC] = MetodoPotencias(A, b, tol, maxiter)
% [l, z, incr, iter, ACOC] = MetodoPotencias(A, b, tol, maxiter)
% Aproxima el valor propio dominante de A
% PARAMETROS:
% A -> Matriz A de la que se quiere sacar el valor propio
% b -> Estimacion del vector propio asociado
% tol -> tolerancia. Detiene ejecucion cuando se alcanza un incremento
% menor
% maxiter -> numero de iteraciones tras las que cesa la ejecucion

    b = b(:); z = b / (b' * b);
    iter = 0; incr = zeros(maxiter + 1, 1); incr(1) = 1;
    ACOC = [];
    p = A * z;

    while iter < maxiter && incr(iter + 1) > tol
        l = z' * p;

        z = p / norm(p);
        p = A * z;
        iter = iter + 1;
        incr(iter + 1) = norm(z - (p / norm(p)));
    end

    if iter >= maxiter
        l = 'No ha convergido';
    else
        incr = incr(2:iter + 1);
        ACOC = fACOC(incr);
    end
end