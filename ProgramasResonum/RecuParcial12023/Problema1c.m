function [sol, incr, iter, ACOC] = Problema1c(A, b, x0, maxiter)
    
    tol = 1e-10;
    w = 1.75;

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
        incr(iter + 1) = norm(A * x1 - b, 2);
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