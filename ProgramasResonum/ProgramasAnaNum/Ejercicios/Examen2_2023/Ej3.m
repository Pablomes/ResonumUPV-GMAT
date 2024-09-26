function [p, c, err, N] = Ej3(f, n, tol)
% [p, c] = AproximacionChebyshevPGorro(f, n)
% Aproxima la funcion a base de los polinomios mejorados P~ de Legendre
% PARAMETROS:
% f -> expresión simbolica a aproximar
% n -> grado
err = [tol + 1]; N = 0;

while err(end) >= tol

    c = sym(zeros(N+1, 1));
    polinomios = sym(zeros(N+1, 1)); syms x;
    for i = 0:N
    polinomios(i+1) = LegendreP(i, x);
    end

    p = 0;

    for i = 0:N
        cK = ((2 * i + 1) / 2) * int(f * polinomios(i+1), x, -1, 1);
    
        c(i+1) = cK;

        p = p + (cK * polinomios(i+1)); 
    end
       
    xe = linspace(-1, 1, n + 1);
    err = [err , max(abs(subs(f, x, xe) - subs(p, x, xe)))];

    N = N + 1;
end

N = N - 1;

end