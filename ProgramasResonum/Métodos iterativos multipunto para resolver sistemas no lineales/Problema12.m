function [sol, incr, incr2, iter, ACOC] = Problema12(x0, n, tol, maxiter)

    x0 = x0(:);
    iter = 0;
    incr = zeros(maxiter + 1, 1); incr(1) = 1;
    incr2 = zeros(maxiter + 1, 1); incr2(1) = 1;
    ACOC = []; [Fx, dFx] = func(x0, n);

    while incr(iter + 1) + incr2(iter + 1) > tol && iter < maxiter
        z = dFx \ Fx;

        y = x0 - z;

        [Fy, dFy] = func(y, n);

        z = dFx \ Fy;

        c = (3 * dFx + dFy) * z;

        z = (5 * dFy - dFx) \ c;

        x1 = y - z;

        iter = iter + 1;
        incr(iter + 1) = norm(x1 - x0);
        [Fx, dFx] = func(x1, n);
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

function [F, dF] = func(x, n)

    F = zeros(n, 1);
    dF = zeros(n, n);

    sumSq = sum(x.^2);

    for i = 1:n
        F(i) = atan(x(i)) + 1 - 2 * (sumSq - n * x(i)^2); 

        for j = 1:n
            if i == j
                dF(i, j) = 1 / (1 + x(j)^2) + 4 * (n - 1) * x(j);
            else
                dF(i, j) = -4 * x(j);
            end
        end
    end
end