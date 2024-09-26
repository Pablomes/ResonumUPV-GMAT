function p = ChebyshevT(n, var)
% p = ChebyshevT(n, var)
% Calcula el polinomio T~ de Chebyshev de grado dado
% PARAMETROS:
% n -> grado
% var -> expresion simbolica (suele ser x)

if n == 0
    p = 1;
elseif n == 1
    p = var;
else

p = expand(2 * var * ChebyshevT(n-1, var) - ChebyshevT(n-2, var));

end
end