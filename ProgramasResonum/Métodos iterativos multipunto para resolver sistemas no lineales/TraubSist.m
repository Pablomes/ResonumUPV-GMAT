function [sol, incr, incr2, iter, ACOC] = TraubSist(F, dF, x0, tol, maxiter)
% [sol, incr, incr2, iter, ACOC] = TraubSist(F, dF, x0, tol, maxiter)
% Aproxima la solucion F(x) = 0 usando el metodo multipaso de Traub
% PARAMETROS:
% F -> Sistema tal que F(x) = 0; Debe ser una funcion anonima.
% dF -> Jacobiana de F(x). Debe ser una matriz anonima.
% x0 -> punto de inicio
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
    z = zeros(n, 1);

    while incr(iter + 1) + incr2(iter + 1) > tol && iter < maxiter
        LU = lu(valdF);
        % SUSTITUCION DIRECTA
        temp(1) = valF(1);

        for i = 2:n
            temp(i) = valF(i) - LU(i, 1:i - 1) * temp(1:i - 1);
        end
        
        % SUSTITUCION INVERSA
        z(n) = temp(n) / LU(n, n);

        for i = (n - 1):-1:1
            z(i) = (temp(i) - LU(i, i + 1:n) * z(i + 1:n)) / LU(i, i);
        end

        y = x0 - z;

        valFy = F(y);

        % SUSTITUCION DIRECTA
        temp(1) = valFy(1);

        for i = 2:n
            temp(i) = valFy(i) - LU(i, 1:i - 1) * temp(1:i - 1);
        end
        
        % SUSTITUCION INVERSA
        z(n) = temp(n) / LU(n, n);

        for i = (n - 1):-1:1
            z(i) = (temp(i) - LU(i, i + 1:n) * z(i + 1:n)) / LU(i, i);
        end

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