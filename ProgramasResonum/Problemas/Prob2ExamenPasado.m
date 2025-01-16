function [sol, incr, incr2, iter, ACOC] = Prob2ExamenPasado(F, dF, x0, maxiter)

    x0 = x0(:); tol = 1e-10;
    iter = 0;
    incr = zeros(maxiter + 1, 1); incr(1) = 1;
    incr2 = zeros(maxiter + 1, 1); incr2(1) = 1;
    ACOC = []; F = @(x) reshape(F(x), [], 1);
    Fx = F(x0); dFx = dF(x0);

    while incr(iter + 1) + incr2(iter + 1) > tol && iter < maxiter
        z = (dFx \ Fx);

        y = x0 - z;

        dFy = dF(y);

        c = dFy \ Fx;

        w = (1/2) * (z + c);

        x1 = x0 - w;

        iter = iter + 1;
        incr(iter + 1) = norm(x1 - x0);
        Fx = F(x1); dFx = dF(x1);
        incr2(iter + 1) = norm(Fx);
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