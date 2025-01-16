function [sol, incr, incr2, iter, ACOC] = Prob16(x0, tol, maxiter)

    x0 = x0(:);
    iter = 0;
    incr = zeros(maxiter + 1, 1); incr(1) = 1;
    incr2 = zeros(maxiter + 1, 1); incr2(1) = 1;
    ACOC = []; [Fx, dFx] = func16(x0);

    while incr(iter + 1) + incr2(iter + 1) > tol && iter < maxiter
        z = dFx \ Fx;
        x1 = x0 - z;
        
        iter = iter + 1;
        incr(iter + 1) = norm(x1 - x0);
        [Fx, dFx] = func16(x1);
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

function [Fx, dFx] = func16(u)
    Fx = zeros(9, 1);
    dFx = zeros(9, 9);

    xk = [0, 0.5, 1, 0, 0.5, 1, 0, 0.5, 1];
    yk = [0, 0, 0, 1, 1, 1, 2, 2, 2];
    pesos = [0, 0, 0, 2, 4, 2, 2, 4, 2];

    for i = 1:9
        Fx(i) = sin(xk(i) * yk(i)) + (exp(xk(i) * yk(i)) / 4) * (u(4).^2 + 2 * (u(5).^2) + u(6).^2 + u(7).^2 + 2 * (u(8).^2) + u(9).^2) - u(i);

        for j = 1:9
            if i == j
                dFx(i, j) = (exp(xk(i) * yk(i)) / 4) * (pesos(j) * u(j)) - 1;
            else
                dFx(i, j) = (exp(xk(i) * yk(i)) / 4) * (pesos(j) * u(j));
            end
        end
    end
end