function [sol, incr, incr2, iter, ACOC] = TrapeciosSist(F, dF, x0, tol, maxiter)
% [sol, incr, incr2, iter, ACOC] = TrapeciosSist(F, dF, x0, tol, maxiter)
% Aproxima la solucion F(x) = 0 usando el metodo multipaso de Trapecios
% PARAMETROS:
% F -> Sistema tal que F(x) = 0; Debe ser una funcion anonima.
% dF -> Jacobiana de F(x). Debe ser una matriz anonima.
% x0 -> punto de inicio
% tol -> tolerancia. Detiene ejecucion cuando se alcanza un incremento
% menor
% maxiter -> numero de iteraciones tras las que cesa la ejecucion   

    x0 = x0(:);
    iter = 0;
    incr = zeros(maxiter + 1, 1); incr(1) = 1;
    incr2 = zeros(maxiter + 1, 1); incr2(1) = 1;
    ACOC = []; F = @(x) reshape(F(x), [], 1);
    valF = F(x0); valdF = dF(x0);
    
    while incr(iter + 1) + incr2(iter + 1) > tol && iter < maxiter
        z = valdF \ valF;

        y = x0 - z;

        z = (valdF + dF(y)) \ valF;

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