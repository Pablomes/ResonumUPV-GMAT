function [p, c] = AproximacionLegendrePGorro(f, n)
% [p, c] = AproximacionChebyshevPGorro(f, n)
% Aproxima la funcion a base de los polinomios mejorados P~ de Legendre
% PARAMETROS:
% f -> expresión simbolica a aproximar
% n -> grado

c = sym(zeros(n+1, 1));
polinomios = sym(zeros(n+1, 1)); syms x;
for i = 0:n
    polinomios(i+1) = LegendreP(i, x);
end

p = 0;

for i = 0:n
    cK = ((2 * i + 1) / 2) * int(f * polinomios(i+1), x, -1, 1);
    
    c(i+1) = cK;

    p = p + (cK * polinomios(i+1)); 
end

fplot(f, [-1, 1])
hold on
fplot(p, [-1, 1])
legend("f(x)", "p" + n + "(x)")
end