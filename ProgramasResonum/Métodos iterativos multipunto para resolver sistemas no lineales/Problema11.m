function [sol, incr, incr2, iter, ACOC] = Problema11(x0, maxiter)

    x0 = x0(:);
    iter = 0; tol = 1e-10;
    incr = zeros(maxiter + 1, 1); incr(1) = 1;
    incr2 = zeros(maxiter + 1, 1); incr2(1) = 1;
    ACOC = []; [Fx, dFx] = func(x0);

    while incr(iter + 1) + incr2(iter + 1) > tol && iter < maxiter
        z = dFx \ Fx;

        y = x0 - z;

        [Fy, dFy] = func(y);

        z = dFy \ (Fx - Fy);

        x1 = x0 - z;

        iter = iter + 1;
        incr(iter + 1) = norm(x1 - x0);
        [Fx, dFx] = func(x1);
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

function [F, dF] = func(x)
    
    F = [x(2) + x(3) + x(4) - x(1) * exp(-x(1)); x(1) + x(3) + x(4) - x(2) * exp(-x(2)); x(1) + x(2) + x(4) - x(3) * exp(-x(3)); x(1) + x(2) + x(3) - x(4) * exp(-x(4))];
    dF = ones(4, 4);

    for i = 1:4
        dF(i, i) = -exp(-x(i)) + x(1) * exp(-x(1));
    end
end