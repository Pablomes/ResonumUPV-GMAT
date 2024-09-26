function p = AproximacionLegendreBasica(f, n)
% p = AproximacionLegendreBasica(f, n)
% Aproxima la funcion a base de los polinomios basicos de Legendre
% PARAMETROS:
% f -> expresión simbolica a aproximar
% n -> grado


p = AproximacionPolinomiosOrtogonales(f, n, 1, -1, 1);

end