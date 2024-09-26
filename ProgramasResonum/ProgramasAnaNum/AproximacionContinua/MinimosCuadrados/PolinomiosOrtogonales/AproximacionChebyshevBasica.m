function p = AproximacionChebyshevBasica(f, n)
% p = AproximacionChebyshevBasica(f, n)
% Aproxima la funcion a base de los polinomios basicos de Chebyshev
% PARAMETROS:
% f -> expresión simbolica a aproximar
% n -> grado

syms x;
p = AproximacionPolinomiosOrtogonales(f, n, 1/sqrt(1 - x^2), -1, 1);

end