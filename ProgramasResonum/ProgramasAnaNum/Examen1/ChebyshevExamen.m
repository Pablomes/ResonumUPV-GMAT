function [p, a, T] = ChebyshevExamen(f, n)
% p = AproximacionChebyshevTGorro(f, n)
% Aproxima la funcion a base de los polinomios mejorados T~ de Chebyshev
% PARAMETROS:
% f -> expresión simbolica a aproximar
% n -> grado

a = sym(zeros(n+1, 1));

polinomios = sym(zeros(n+1, 1)); syms x;
for i = 0:n
    polinomios(i+1) = ChebyshevT(i, x);
end

p = 0;

a0 = (2/pi) * int((f * polinomios(1)) / (sqrt(1 - x^2)), x, -1, 1);
a(1) = a0;
p = p + (a0 / 2);

for i = 1:n
    a0 = (2/pi) * int((f * polinomios(i+1)) / (sqrt(1 - x^2)), x, -1, 1);
    a(i + 1) = a0;

    p = p + a0 * polinomios(i+1); 
end

T = polinomios;

fplot(f, [-1, 1])
hold on
fplot(p, [-1, 1])
legend("f(x)", "p" + n + "(x)")
end