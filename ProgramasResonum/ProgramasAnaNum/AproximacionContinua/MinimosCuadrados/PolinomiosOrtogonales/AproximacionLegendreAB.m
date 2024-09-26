function [p, c] = AproximacionLegendreAB(f, n, a, b)
% [p, c] = AproximacionLegendreAB(f, n, a, b)
% Aproxima la funcion a base de los polinomios mejorados P~ de Legendre en
% el intervalo [a, b]
% PARAMETROS:
% f -> expresión simbolica a aproximar
% n -> grado
% [a, b] -> Intervalo de aproximacion

c = sym(zeros(n+1, 1));

polinomios = sym(zeros(n+1, 1)); syms x; syms t;

ft = subs(f, x, ((b-a) * t + (b + a)) / 2);

for i = 0:n
    polinomios(i+1) = LegendreP(i, t);
end

p = 0;

for i = 0:n
    cK = ((2 * i + 1) / 2) * int(ft * polinomios(i+1), t, -1, 1);
    c(i+1) = subs(cK, t, (2 * x - b - a) / (b - a));

    p = p + (cK * polinomios(i+1)); 
end

p = subs(p, t, (2 * x - b - a) / (b - a));

fplot(f, [a, b])
hold on
fplot(p, [a, b])
legend("f(x)", "p" + n + "(x)")
end