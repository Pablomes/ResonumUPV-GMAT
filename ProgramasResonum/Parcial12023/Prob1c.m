function [sol, incr, iter, ACOC] = Prob1c(A, b, x0, maxiter)   

    tol = 1e-10;

    b = b(:); x0 = x0(:); n = size(A, 1);
    incr = zeros(maxiter + 1, 1); incr(1) = 1;
    iter = 0;
    ACOC = [];

    while iter < maxiter && incr(iter + 1) > tol
        r = b - A * x0;
        t = (r' * r) / (r' * A * r);

        x1 = x0 + t * r;
        
        iter = iter + 1;
        incr(iter + 1) = sqrt(r' * r);
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