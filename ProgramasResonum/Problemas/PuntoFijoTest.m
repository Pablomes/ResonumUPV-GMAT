function [sol, incr, iter, ACOC] = PuntoFijoTest(x0, tol, maxiter)
    
    iter = 0; ACOC = [];
    incr = zeros(maxiter, 1);
    incr(1) = 1;

    G = @(x) [5 / (x(1) - x(2)); 1 - 3 * x(1)];

    while incr(iter + 1) > tol && iter < maxiter
        x1 = G(x0);

        iter = iter + 1;
        incr(iter + 1) = norm(x1 - x0);
        x0 = x1;
    end

    if iter >= maxiter
        sol = 'NO';
    else
        sol = x0;
        incr = incr(2:iter + 1);
        ACOC = fACOC(incr);
end