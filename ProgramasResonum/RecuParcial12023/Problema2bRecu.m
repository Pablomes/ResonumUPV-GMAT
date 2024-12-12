function [R, valP, incr, iter, ACOC] = Problema2bRecu(A, tol, maxiter)

    [~, R] = TridiagonalizacionHouseholder(A);
    iter = 0; incr = zeros(maxiter + 1, 1); incr(1) = 1;
    ACOC = []; n = size(A, 1);

    while iter < maxiter && incr(iter + 1) > tol
        [Q, R2] = qr(R);
        R = R2 * Q;

        iter = iter + 1;
        incr(iter + 1) = norm(R - diag(diag(R)));
    end

    if iter >= maxiter
        valP = "No ha convergido";
    else
        valP = diag(R);
        incr = incr(2:iter + 1);
        ACOC = fACOC(incr);
    end
end