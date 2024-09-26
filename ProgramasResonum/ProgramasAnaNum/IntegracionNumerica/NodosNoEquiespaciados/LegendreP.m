function p = LegendreP(n, var)
% p = LegendreP(n, var)
% Calcula el polinomio P~ de Legendre de grado dado
% PARAMETROS:
% n -> grado
% var -> expresion simbolica (suele ser x)

if n == 0
    p = 1;
elseif n == 1
    p = var;
else

k = n - 1;

p = expand(((2 * k + 1) / (k + 1)) * var * LegendreP(k, var) - ((k) / (k + 1)) * LegendreP(k - 1, var));
end
end