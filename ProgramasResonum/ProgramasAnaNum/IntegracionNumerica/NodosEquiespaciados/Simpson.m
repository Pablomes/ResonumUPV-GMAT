function S = Simpson(f, a, b, n)
% S = Simpson(f, a, b, n)
% Aproxima la integral a base del metodo de Simpson
% PARAMETROS:
% f -> funcion anonima a integrar
% [a, b] -> dominio
% n -> numero de curvas (debe ser par)

if mod(n, 2) ~= 0
    error("n debe ser par!")
end

S = 0;
h = (b - a) / n;

for i = 0:2:(n - 2)
    S = S + f(a + h * i) + 4 * f(a + h * (i + 1)) + f(a + h * (i + 2));
end

S = S * (h / 3);

end