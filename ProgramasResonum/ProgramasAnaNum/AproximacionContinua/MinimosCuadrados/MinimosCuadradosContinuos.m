function [p, coef] = MinimosCuadradosContinuos(f, n, a, b)
% p = MinimosCuadradosContinuos(f, n, a, b)
% Aproxima la funcion a base del metodo de minimos cuadrados
% PARAMETROS:
% f -> expresión simbolica a aproximar
% n -> grado
% [a, b] -> dominio

A = zeros(n+1, n+1); s = zeros(n+1, 1); syms x;

for j = 0:n
    for k = 0:n
        A(j+1, k+1) = (b^(j + k + 1) - a^(j + k + 1)) / (j + k + 1);
    end

    s(j + 1) = int((x^j) * f, x, a, b);
end

coef = Cramer(A, s);

p = 0;

for i = 1:length(coef)
    p = p + coef(i) * x^(i-1);
end

fplot(f, [a, b])
hold on
fplot(p, [a, b])
legend("f(x)", "p" + n + "(x)")

coef = coef(end:-1:1);
coef = coef(:);

end