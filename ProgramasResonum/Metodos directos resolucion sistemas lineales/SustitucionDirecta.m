function sol = SustitucionDirecta(L, b)
% sol = SustitucionDirecta(L, b)
% Devuelve la solución al sistema usando sustitucion directa
% PARAMETROS:
% L -> Matriz del sistema, debe ser cuadrada y triangular inferior
% b -> vector independiente del sistema
    
    n = size(L, 2);
    b = b(:);
    sol = zeros(n, 1);

    sol(1) = b(1);

    for i = 2:n
        sol(i) = b(i) - L(i, 1:i - 1) * sol(1:i - 1);
    end
end