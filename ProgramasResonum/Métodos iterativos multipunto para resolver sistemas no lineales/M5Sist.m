function [sol, incr, incr2, iter, ACOC] = M5Sist(F, dF, x0, alpha, beta, gamma, rho, mu, sigma, tol, maxiter)
% [sol, incr, incr2, iter, ACOC] = M5Sist(F, dF, x0, alpha, beta, gamma, rho, mu, sigma, tol, maxiter)
% Aproxima la solucion F(x) = 0 usando el metodo multipaso de la familia
% HMT
% Nota: El método tiene orden de convergencia 5 cuando alpha = 1,
% beta = rho = 0, gamma = 5, mu = -16/5, sigma = 1/5
% PARAMETROS:
% F -> Sistema tal que F(x) = 0; Debe ser una funcion anonima.
% dF -> Jacobiana de F(x). Debe ser una matriz anonima.
% x0 -> punto de inicio
% alpha, beta, gamma, rho, mu, sigma -> parametros del método
% tol -> tolerancia. Detiene ejecucion cuando se alcanza un incremento
% menor
% maxiter -> numero de iteraciones tras las que cesa la ejecucion  

    x0 = x0(:); [n, ~] = size(x0);
    iter = 0;
    incr = zeros(maxiter + 1, 1); incr(1) = 1;
    incr2 = zeros(maxiter + 1, 1); incr2(1) = 1;
    ACOC = []; F = @(x) reshape(F(x), [], 1);
    valF = F(x0); valdF = dF(x0);
    temp = zeros(n, 1);
    s = zeros(n, 1);

    while incr(iter + 1) + incr2(iter + 1) > tol && iter < maxiter
        LU = lu(valdF);

        % yk = xk - ...
        % SUSTITUCION DIRECTA
        temp(1) = valF(1);

        for i = 2:n
            temp(i) = valF(i) - LU(i, 1:i - 1) * temp(1:i - 1);
        end
        
        % SUSTITUCION INVERSA
        s(n) = temp(n) / LU(n, n);

        for i = (n - 1):-1:1
            s(i) = (temp(i) - LU(i, i + 1:n) * s(i + 1:n)) / LU(i, i);
        end

        y = x0 - alpha * s;

        valFy = F(y);

        %zk = yk - ...
        temp = zeros(n, 1);
        s = zeros(n, 1);
        % SUSTITUCION DIRECTA
        b = (beta * valF + gamma * valFy);
        temp(1) = b(1);

        for i = 2:n
            temp(i) = b(i) - LU(i, 1:i - 1) * temp(1:i - 1);
        end
        
        % SUSTITUCION INVERSA
        s(n) = temp(n) / LU(n, n);

        for i = (n - 1):-1:1
            s(i) = (temp(i) - LU(i, i + 1:n) * s(i + 1:n)) / LU(i, i);
        end

        z = y - s;

        valFz = F(z);

        % x(k+1) = zk - ...
        temp = zeros(n, 1);
        s = zeros(n, 1);
        % SUSTITUCION DIRECTA
        b = (rho * valF + mu * valFy + sigma * valFz);
        temp(1) = b(1);

        for i = 2:n
            temp(i) = b(i) - LU(i, 1:i - 1) * temp(1:i - 1);
        end
        
        % SUSTITUCION INVERSA
        s(n) = temp(n) / LU(n, n);

        for i = (n - 1):-1:1
            s(i) = (temp(i) - LU(i, i + 1:n) * s(i + 1:n)) / LU(i, i);
        end

        x1 = z - s;

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