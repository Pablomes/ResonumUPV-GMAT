function [sol, incr, incr2, iter, ACOC] = Problema9Trapecios(f, fy, fdy, N, tol, maxiter)
    
    a = 0; b = pi; alpha = 1; beta = 2;
    h = (b - a) / (N - 1); yk = (alpha + (1:N) * (beta - alpha) * h / (b - a))';
    iter = 0; ACOC = [];
    incr = zeros(maxiter + 1, 1); incr(1) = 1;
    incr2 = zeros(maxiter + 1, 1); incr2(1) = 1;
    mDiag = zeros(N, 1); zk = zeros(N, 1);
    lDiag = zeros(N - 1, 1); uDiag = zeros(N - 1, 1);
    bk = zeros(N, 1);

    xk = (a : h : b)';

    zk(1) = (yk(2) - h * yk(1)) / h;
    mDiag(1) = (2 * h - 2) + (h^2) * (fy(xk(1), yk(1), zk(1)) - fdy(xk(1), yk(1), zk(1)));
    uDiag(1) = 2 + h * fdy(xk(1), yk(1), zk(1));
    bk(1) = (2 * h - 2) * yk(1) + 2 * yk(2) + h^2 * f(xk(1), yk(1), zk(1)) - 2 * h;

    for i = 2 : (N - 1)
        zk(i) = (yk(i + 1) - yk(i - 1)) / (2 * h);
        mDiag(i) = -2 + h^2 * fy(xk(i), yk(i), zk(i));
        uDiag(i) = 1 + (h/2) * fdy(xk(i), yk(i), zk(i));
        lDiag(i - 1) = 1 - (h/2) * fdy(xk(i), yk(i), zk(i));
        bk(i) = yk(i - 1) - 2 * yk(i) + yk(i + 1) + h^2 * f(xk(i), yk(i), zk(i));
    end

    zk(N) = (yk(N) - 2) / 2;
    mDiag(N) = (h - 2) + h^2 * (fy(xk(N), yk(N), zk(N)) + (1/2) * fdy(xk(N), yk(N), zk(N)));
    lDiag(N - 1) = 2;
    bk(N) = 2 * yk(N - 1) + (h - 2) * yk(N) + h^2 * f(xk(N), yk(N), zk(N)) - 2 * h;

    while incr(iter + 1) + incr2(iter + 1) > tol && iter < maxiter
        z = ResolucionDirectaCrout(mDiag, lDiag, uDiag, bk);

        y = yk - z;
        
        zk(1) = (y(2) - h * y(1)) / h;
        mDiag(1) = mDiag(1) + ((2 * h - 2) + (h^2) * (fy(xk(1), y(1), zk(1)) - fdy(xk(1), y(1), zk(1))));
        uDiag(1) = uDiag(1) + 2 + h * fdy(xk(1), y(1), zk(1));

        for i = 2 : (N - 1)
            zk(i) = (y(i + 1) - y(i - 1)) / (2 * h);
            mDiag(i) = mDiag(i) -2 + h^2 * fy(xk(i), y(i), zk(i));
            uDiag(i) = uDiag(i) + 1 + (h/2) * fdy(xk(i), y(i), zk(i));
            lDiag(i - 1) = lDiag(i - 1) + 1 - (h/2) * fdy(xk(i), y(i), zk(i));
        end

        zk(N) = (y(N) - 2) / 2;
        mDiag(N) = mDiag(N) + (h - 2) + h^2 * (fy(xk(N), y(N), zk(N)) + (1/2) * fdy(xk(N), y(N), zk(N)));
        lDiag(N - 1) = lDiag(N - 1) + 2;

        z = ResolucionDirectaCrout(mDiag, lDiag, uDiag, bk);

        y1 = yk - z;

        iter = iter + 1;
        incr(iter + 1) = norm(y1 - yk);

        yk = y1;

        zk(1) = (yk(2) - h * yk(1)) / h;
        mDiag(1) = (2 * h - 2) + (h^2) * (fy(xk(1), yk(1), zk(1)) - fdy(xk(1), yk(1), zk(1)));
        uDiag(1) = 2 + h * fdy(xk(1), yk(1), zk(1));
        bk(1) = (2 * h - 2) * yk(1) + 2 * yk(2) + h^2 * f(xk(1), yk(1), zk(1)) - 2 * h;

        for i = 2 : (N - 1)
            zk(i) = (yk(i + 1) - yk(i - 1)) / (2 * h);
            mDiag(i) = -2 + h^2 * fy(xk(i), yk(i), zk(i));
            uDiag(i) = 1 + (h/2) * fdy(xk(i), yk(i), zk(i));
            lDiag(i - 1) = 1 - (h/2) * fdy(xk(i), yk(i), zk(i));
            bk(i) = yk(i - 1) - 2 * yk(i) + yk(i + 1) + h^2 * f(xk(i), yk(i), zk(i));
        end

        zk(N) = (yk(N) - 2) / 2;
        mDiag(N) = (h - 2) + h^2 * (fy(xk(N), yk(N), zk(N)) + (1/2) * fdy(xk(N), yk(N), zk(N)));
        lDiag(N - 1) = 2;
        bk(N) = 2 * y(N - 1) + (h - 2) * y(N) + h^2 * f(xk(N), yk(N), zk(N)) - 2 * h;

        incr2(iter + 1) = norm(bk);
    end

    if iter >= maxiter
        sol = 'No ha convergido';
    else
        sol = yk;
        incr = incr(2:iter + 1);
        ACOC = fACOC(incr);
    end

end