function [S, coef] = SplineCubicoFrontera(xi, fi, id, fd)
% [S, coef] = SplineCubicoFrontera(xi, fi, id, fd)
% Aproxima los puntos usando un spline con las primeras derivadas de los
% puntos extremos dada
% PARAMETROS:
% xi -> lista de coordenadas x
% fi -> lista de coordenadas y
% id -> primera derivada del extremo izquierdo del spline
% fd -> primera derivada del extremo derecho del spline

xi = xi(:); fi = fi(:); n = length(xi); 
S = sym(zeros(n-1, 1)); syms x;

border = ((xi(n) - xi(1)) / n) / 4;
scatter(xi, fi)
hold on

A = [];
b = [];

% AÑADE LAS ECUACIONES PARA QUE LOS SPLINES TOMEN LOS VALORES DE LOS PUNTO
for i = 1:(n-1)
    A = [A;
        [zeros(1, 4 * (i - 1)),[1, 0, 0, 0], zeros(1, 4*(n-1-i))]; 
        [zeros(1, 4 * (i - 1)),[1, (xi(i+1)-xi(i)), (xi(i+1)-xi(i))*(xi(i+1)-xi(i)), (xi(i+1)-xi(i))*(xi(i+1)-xi(i))*(xi(i+1)-xi(i))], zeros(1, 4*(n-1-i))]];
    b = [b;
        fi(i);
        fi(i+1)];
end
% AÑADE LAS ECUACIONES PARA QUE LOS SPLINES TOMEN EL MISMO VALOR DE PRIMERA
% Y SEGUNDA DERIVADA EN LOS PUNTOS
for i = 1:(n-2)
    A = [A;
        [zeros(1, 4 * (i - 1)), [0, 1, 2 * (xi(i+1)-xi(i)), 3 * (xi(i+1)-xi(i)) * (xi(i+1)-xi(i))], [0, -1, 0, 0], zeros(1, 4 * (n-2-i))];
        [zeros(1, 4 * (i - 1)), [0, 0, 2, 6 * (xi(i+1)-xi(i))], [0, 0, -2, 0], zeros(1, 4 * (n-2-i))]];
    b = [b;
        0;
        0];
end

% AÑADE LAS CONDICIONES INICIALES PARA LA PRIMERA DERIVADA

A = [A;
    [[0, 1, 0, 0], zeros(1, 4 * (n-2))];
    [zeros(1, 4 * (n-2)), [0, 1, 2 * (xi(n)-xi(n-1)), 6 * (xi(n)-xi(n-1)) * (xi(n)-xi(n-1))]]];
b = [b;
    id;
    fd];

% RESOLVEMOS EL SISTEMA

c = Cramer(A, b);

% GENERAMOS LOS SPLINES USANDO LOS COEFICIENTES CALCULADOS
coef = [];

for i = 1:(n-1)
    coef = [coef; [c(4 * (i-1) + 4), c(4 * (i-1) + 3), c(4 * (i-1) + 2), c(4 * (i-1) + 1)]];
    S(i) = c(4 * (i-1) + 1) + c(4 * (i-1) + 2) * (x - xi(i)) + c(4 * (i-1) + 3) * (x - xi(i)) * (x - xi(i)) + c(4 * (i-1) + 4) * (x - xi(i)) * (x - xi(i)) * (x - xi(i));
    fplot(S(i), [xi(i), xi(i+1)])
    hold on
end
end