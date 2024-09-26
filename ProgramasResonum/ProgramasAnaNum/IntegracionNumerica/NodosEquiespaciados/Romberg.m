function R = Romberg(f, a, b, n, dim)
% R = Romberg(f, a, b, n), dim
% Aproxima la integral usando la matriz de Romberg y el metodo de trapecios
% PARAMETROS:
% f -> funcion anonima a integrar
% [a, b] -> dominio
% n -> numero de subintervalos iniciales
% dim -> dimension de la matriz cuadrada de Romberg

R = zeros(dim, dim); h = (b - a) / n;

for i = 0:dim - 1
    R(i + 1, 1) = Trapecios(f, a, b, n * 2^(i));
end

for i = 1:dim - 1
    for j = 1:dim - i
        R(j, i + 1) = (1 / (4^(i) - 1)) * (4^(i) * R(j + 1, i) - R(j, i));
    end
end
end