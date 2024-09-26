function p = AproximacionChebyshevAB(f, n, a, b)
% p = AproximacionChebyshevAB(f, n, a, b)
% Aproxima la funcion a base de los polinomios mejorados T~ de Chebyshev en
% el intervalo [a, b]
% PARAMETROS:
% f -> expresión simbolica a aproximar
% n -> grado
% [a, b] -> Intervalo de aproximacion

polinomios = sym(zeros(n+1, 1)); syms x t;

ft = subs(f, x, ((b-a) * t + (b + a)) / 2);

for i = 0:n
    polinomios(i+1) = ChebyshevT(i, t)
end

p = 0;

a0 = (2/pi) * int((ft * polinomios(1)) / (sqrt(1 - t^2)), t, -1, 1);
p = p + (a0 / 2);

for i = 1:n
    a0 = (2/pi) * int((ft * polinomios(i+1)) / (sqrt(1 - t^2)), t, -1, 1);

    p = p + a0 * polinomios(i+1); 
end

p = subs(p, t, (2 * x - b - a) / (b - a));

fplot(f, [a, b])
hold on
fplot(p, [a, b])
legend("f(x)", "p" + n + "(x)")
end