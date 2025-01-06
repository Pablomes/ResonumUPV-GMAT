function [sol, incr, iter, ACOC] = DifFinNoLinNewton(f, fy, fdy, a, b, alpha, beta, N, tol, maxiter)
% [sol, incr, iter, ACOC] = DifFinNoLinNewton(f, fy, fdy, a, b, alpha, beta, tol, maxiter)
% Aproxima la solucion F(x) = 0 usando el metodo multipaso de Jarratt
% PARAMETROS:
% f -> función tal que y'' = f(x, y, y'). Debe ser funcion anonima f = @(x, y, z)
% fy -> derivada de f respecto a y. Debe ser funcion anonima fy = @(x, y, z)
% fdy -> derivada de f respecto de y'. Debe ser funcion anonima fdy = @(x, y, z)
% [a, b] -> Intervalo del problema
% alpha -> condicion inicial tal que f(a) = alpha
% beta -> condicion inicial tal que f(b) = beta
% N -> numero de incognitas
% tol -> tolerancia. Detiene ejecucion cuando se alcanza un incremento
% menor
% maxiter -> numero de iteraciones tras las que cesa la ejecucion  

    h = (b - a) / (N + 1); yk = (alpha + (1:N) * (beta - alpha) * h / (b - a))';
    iter = 0; ACOC = [];
    incr = zeros(maxiter + 1, 1); incr(1) = 1;
    diag = zeros(N, 1); xk = zeros(N, 1); zk = zeros(N, 1);
    lDiag = zeros(N - 1, 1); uDiag = zeros(N - 1, 1);
    bk = zeros(N, 1);

    while incr(iter + 1) > tol && iter < maxiter

        xk(1) = a + h; zk(1) = (yk(2) - alpha) / (2 * h);
        diag(1) = 2 + h^2 * fy(xk(1), yk(1), zk(1));
        uDiag(1) = -1 + (h/2) * fdy(xk(1), yk(1), zk(1));
        bk(1) = -(2 * yk(1) - yk(2) + h^2 * f(xk(1), yk(1), zk(1)) - alpha);

        for i = 2 : (N - 1)
            xk(i) = a + i * h; zk(i) = (yk(i + 1) - yk(i - 1)) / (2 * h);

            diag(i) = 2 + h^2 * fy(xk(i), yk(i), zk(i));
            uDiag(i) = -1 + (h/2) * fdy(xk(i), yk(i), zk(i));
            lDiag(i - 1) = -1 - (h/2) * fdy(xk(i), yk(i), zk(i));
            bk(i) = -(2 * yk(i) - yk(i + 1) - yk(i - 1) + h^2 * f(xk(i), yk(i), zk(i)));
        end


        xk(N) = b - h; zk(N) = (beta - yk(N - 1)) / (2 * h);
        diag(N) = 2 + h^2 * fy(xk(N), yk(N), zk(N));
        lDiag(N - 1) = -1 - (h/2) * fdy(xk(N), yk(N), zk(N));
        bk(N) = -(2 * yk(N) - yk(N - 1) + h^2 * f(xk(N), yk(N), zk(N)) - beta);

        z = ResolucionDirectaCrout(diag, lDiag, uDiag, bk);

        yk = yk + z;

        iter = iter + 1;
        incr(iter + 1) = norm(z);
    end

    if iter >= maxiter
        sol = 'No ha convergido';
    else
        sol = yk;
        incr = incr(2:iter + 1);
        ACOC = fACOC(incr);
    end

end