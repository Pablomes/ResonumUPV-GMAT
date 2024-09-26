function p = PolinomioHermite(xi, fi, dfi)
% p = PolinomioHermite(xi, fi, dfi)
% Aproxima los puntos calculando el polinomio de Hermite
% PARAMETROS:
% xi -> lista de coordenadas x
% fi -> lista de coordenadas y
% dfi -> lista de primeras derivadas en las coordenadas (x, y)

xi = xi(:); fi = fi(:); dfi = dfi(:); n = length(xi);
syms x;

p = 0;

for i = 1:n
    L = 1;
    for j = 1:n
        if i == j
            L = L * 1;
        else
            L = L * ((x - xi(j)) / (xi(i) - xi(j)));
        end
    end

    dL = diff(L);
    dLEval = subs(dL, x, xi(i));

    H = (1 - 2 * (x - xi(i)) * dLEval) * (L * L);
    Hgorro = (x - xi(i)) * (L * L);

    p = p + ((fi(i) * H) + dfi(i) * Hgorro);  
end

border = ((xi(n) - xi(1)) / n) / 4;
fplot(p, [xi(1) - border, xi(n) + border])
hold on
scatter(xi, fi)

end