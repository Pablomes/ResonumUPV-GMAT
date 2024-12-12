function [valP, incr, iter, ACOC] = QRFrancis(A, tol, maxiter)
% [valP, incr, iter, ACOC] = QRFrancis(A tol, maxiter)
% Aproxima los valores propios de A
% PARAMETROS:
% A -> Matriz A simetrica de la que se quiere sacar los valores propios
% tol -> tolerancia. Detiene ejecucion cuando se alcanza un incremento
% menor
% maxiter -> numero de iteraciones tras las que cesa la ejecucion

    [~, R] = TridiagonalizacionHouseholder(A); % NO CAMBIAR ESTE, CAMBIAR EL OTRO
    iter = 0; incr = zeros(maxiter + 1, 1); incr(1) = 1;
    ACOC = []; n = size(A, 1);

    while iter < maxiter && incr(iter + 1) > tol
        [Q, R2] = qr(R);
        R = R2 * Q;

        sum = 0;

        for i = 1:n
            for j = 1:n
                if i ~= j
                    sum = sum + abs(R(i, j));
                end
            end
        end

        iter = iter + 1;
        incr(iter + 1) = sum;
    end

    if iter >= maxiter
        valP = "No ha convergido";
    else
        valP = diag(R);
        incr = incr(2:iter + 1);
        ACOC = fACOC(incr);
    end
end