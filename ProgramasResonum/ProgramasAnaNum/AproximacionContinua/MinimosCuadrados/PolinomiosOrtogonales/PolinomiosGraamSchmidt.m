function [p] = PolinomiosGraamSchmidt(n, w, a, b)
% p = PolinomiosGraamSchmidt(n, w, a, b)
% Utiliza el Tma de Graam Schmidt para calcular los n polinomios
% ortogonales de la funcion peso dada
% PARAMETROS:
% n -> grado
% w -> funcion peso
% [a, b] -> dominio

p = sym(zeros(n+1, 1)); syms x;
alpha = sym(zeros(n));
p(1) = 1;
alpha(1) = int(w * p(1) * p(1), x, a, b);
p(2) = x - ((1/alpha(1)) * int(x * w * p(1) * p(1), x, a, b));

for i = 2:n
    alpha(i) = int(w * p(i) * p(i), x, a, b);
    bK = ((1/alpha(i)) * int(x * w * p(i) * p(i), x, a, b));
    cK = ((1/alpha(i-1)) * int(x * w * p(i) * p(i-1), x, a, b));
    p(i+1) = expand((x - bK) * p(i) - cK * p(i-1));
end
end