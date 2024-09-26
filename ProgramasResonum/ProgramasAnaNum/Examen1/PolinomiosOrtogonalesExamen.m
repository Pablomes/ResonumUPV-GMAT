function p = PolinomiosOrtogonalesExamen(f, n, w, a, b, polinomios)
% p = AproximacionPolinomiosOrtogonales(f, n, w, a, b)
% Aproxima la funcion a base de la funcion peso dada en el dominio dado
% PARAMETROS:
% f -> expresión simbolica a aproximar
% n -> grado
% w -> funcion peso
% [a, b] -> dominio

%polinomios = PolinomiosGraamSchmidt(n, w, a, b);
p = 0; syms x;

for i = 0:n
    alphaJ = int(w * polinomios(i+1) * polinomios(i+1), x, a, b);
    aJ = ((1/alphaJ) * int(w * f * polinomios(i+1), x, a, b));

    p = p + (aJ * polinomios(i+1)); 
end

fplot(f, [a, b])
hold on
fplot(p, [a, b])
legend("f(x)", "p" + n + "(x)")

end