function [sol, incr, incr2, iter, ACOC] = JarrattSist(F, dF, x0, tol, maxiter)

    x0 = x0(:);
    iter = 0;
    incr = zeros(maxiter + 1, 1); incr(1) = 1;
    incr2 = zeros(maxiter + 1, 1); incr2(1) = 1;
    ACOC = []; F = @(x) reshape(F(x), [], 1);
    valF = F(x0); valdF = dF(x0);
    
    while incr(iter + 1) + incr2(iter + 1) > tol && iter < maxiter
        z = (2/3) * (valdF \ valF);

        y = x0 - z;

        z = valdF \ valF;

        valdFy = 3 * dF(y);

        c = (valdFy + valdF) * z;

        z = (1/2) * ((valdFy - valdF) \ c);

        x1 = y - z;

        iter = iter + 1;
        incr(iter + 1) = norm(x1 - x0);
        valF = F(x1); valdF = dF(x1);
        incr2(iter + 1) = norm(valF);
        x0 = x1;
    end

    if iter >= maxiter
        sol = 'No ha convergido';
    else
        sol = x0;
        incr = incr(2:iter + 1);
        incr2 = incr2(2: iter + 1);
        ACOC = fACOC(incr);
    end

end