function T = Trapecios(f, a, b, n)
% T = Trapecios(f, a, b, n)
% Aproxima la integral a base de trapecios
% PARAMETROS:
% f -> funcion anonima a integrar
% [a, b] -> dominio
% n -> numero de trapecios

T = 0;
h = (b - a) / n;

for i = 0:(n-1)
    T = T + f(a + h * i) + f(a + h * (i + 1));
end

T = T * (h / 2);

end