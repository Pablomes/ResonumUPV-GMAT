function [sol, incr, incr2, iter, ACOC] = Problema12Jarratt(x0, n, tol, maxiter)
% [sol, incr, incr2, iter, ACOC] = JarrattSist(F, dF, x0, tol, maxiter)
% Aproxima la solucion F(x) = 0 usando el metodo multipaso de Jarratt
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
    ACOC = []; 
    [valF, valdF] = func(x0, n);
    
    while incr(iter + 1) + incr2(iter + 1) > tol && iter < maxiter
        z = (2/3) * (valdF \ valF);

        y = x0 - z;

        z = valdF \ valF;
        
        [~, temp] = func(y, n);

        valdFy = 3 * temp;

        c = (valdFy + valdF) * z;

        z = (1/2) * ((valdFy - valdF) \ c);

        x1 = x0 - z;

        iter = iter + 1;
        incr(iter + 1) = norm(x1 - x0);
        [valF, valdF] = func(x1, n);
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