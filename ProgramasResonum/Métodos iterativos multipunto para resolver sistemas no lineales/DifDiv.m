function DD = DifDiv(F, x, y)
% DD = DifDiv(F, x, y)
% Retorna el operador de diferencias divididas [x, y; F]
% PARAMETROS:
% F -> Matriz que define el sistema. Debe ser una funcion anonima.
% x, y -> vectores para la diferencia dividida

    m = length(x);

    DD = zeros(m, m);

    for j = 1 : m
        v1 = x; v1(1:j) = y(1:j);
        v2 = x; v2(1:j - 1) = y(1:j - 1);
        Fv1 = F(v1);
        Fv2 = F(v2);
        DD(:, j) = (Fv1(:) - Fv2(:)) / (y(j) - x(j));
    end
end