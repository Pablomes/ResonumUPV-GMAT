function [p, b, D] = PolinomioNewton(xi, fi)
% [p, b, D] = PolinomioNewton(xi, fi)
% Aproxima los puntos calculando el polinomio de Newton
% PARAMETROS:
% xi -> lista de coordenadas x
% fi -> lista de coordenadas y

xi = xi(:); fi = fi(:); n = length(xi);
syms x; D = sym([fi, zeros(n, n-1)]);

p = fi(1);
b = sym(fi(1));

for col = 2:(n)
    vx = 1;
    for i = 1:col-1
        vx = vx * (x - xi(i));
    end
    for row = col:(n)
        D(row, col) = (D(row, col-1) - D(row-1, col-1))/(xi(row) - xi(row - col + 1));
    end
    b = [b; D(col, col)];
    p = p + (vx * D(col, col));
end

border = ((xi(n) - xi(1)) / n) / 4;
fplot(p, [xi(1) - border, xi(n) + border])
hold on
scatter(xi, fi)

p = expand(p);

%b = b(2:end);

end