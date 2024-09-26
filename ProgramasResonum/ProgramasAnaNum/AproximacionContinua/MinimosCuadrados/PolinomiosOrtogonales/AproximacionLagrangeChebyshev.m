function p = AproximacionLagrangeChebyshev(f, n, a, b)
% p = AproximacionLagrangeChebyshev(f, n, a, b)
% Aproxima la funcion con Lagrange utilizando las raices de Chebyshev como
% nodos
% PARAMETROS:
% f -> expresión simbolica a aproximar
% n -> grado
% [a, b] -> dominio

syms x;

T = ChebyshevT(n+1, x);
roots = solve(T==0);

nodes = zeros(length(roots), 1);
images = zeros(length(nodes), 1);

for i = 1:length(nodes)
    nodes(i) = (b - a) * roots(i) / 2 + (b + a) / 2;
    images(i) = subs(f, x, nodes(i));
end

xi = nodes(:); fi = images(:); n = length(xi);
syms x; L = sym(zeros(n, 1)); p = 0;
for i = 1:n
    e = 1;
    for j = 1:n
        if i == j
            e = e * 1;
        else
            e = e * ((x - xi(j))/(xi(i) - xi(j)));
        end
    end
    L(i) = e;
    p = p + (L(i) * fi(i));
end

fplot(f, [a, b])
hold on
fplot(p, [a, b])
legend("f(x)", "p" + n + "(x)")

end