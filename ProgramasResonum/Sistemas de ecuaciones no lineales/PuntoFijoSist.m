function [sol, incr, iter, ACOC] = PuntoFijoSist(G, x0, tol, maxiter)
% [sol, incr, iter, ACOC] = PuntoFijoSist(G, x0, tol, maxiter)
% Resuelve por el metodo de punto fijo
% PARAMETROS:
% G -> función tal que x = G(x). Debe ser una funcion anonima.
% x0 -> x inicial vector
% tol -> tolerancia
% maxiter -> maximo numero de iteraciones

    x0 = x0(:);
    iter = 0;
    incr = zeros(maxiter + 1, 1); incr(1) = 1;
    ACOC = []; G = @(x) reshape(G(x), [], 1);

    while incr(iter + 1) > tol && iter < maxiter
        x1 = G(x0);

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