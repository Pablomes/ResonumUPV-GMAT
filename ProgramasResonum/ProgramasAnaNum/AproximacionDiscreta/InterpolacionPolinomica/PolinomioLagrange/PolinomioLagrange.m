function p = PolinomioLagrange(xi, fi)
% p = PolinomioLagrange(xi, fi)
% Aproxima los puntos calculando el polinomio de Lagrange
% PARAMETROS:
% xi -> lista de coordenadas x
% fi -> lista de coordenadas y

xi = xi(:); fi = fi(:); n = length(xi);
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

border = ((xi(n) - xi(1)) / n) / 4;
fplot(p, [xi(1) - border, xi(n) + border])
hold on
scatter(xi, fi)

p = expand(p);

end