function p = AproximacionPade(f, n, m, a, b)
% p = AproximacionPade(f, n, m, a, b)
% USA PADE PARA APROXIMAR POR FUNCIONES RACIONALES
% PARAMETROS:
% f -> expresión simbolica a aproximar
% n -> grado de p(x)
% m -> grado de q(x)
% [a, b] -> dominio

syms x;

p = pade(f, x, 0, 'Order', [n, m]);

fplot(f, [a, b])
hold on
fplot(p, [a, b])
legend("f(x)", "p" + n + "(x)")

end