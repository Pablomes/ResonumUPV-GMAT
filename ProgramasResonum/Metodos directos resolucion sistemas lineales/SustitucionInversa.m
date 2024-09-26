function sol = SustitucionInversa(U, z)
% sol = SustitucionInversa(U, z)
% Encuentra la solucion del sistema usando la sustitucion inversa
% PARAMETROS:
% U -> Matriz del sistema, debe ser triangular superior
% z -> Vector independiente del sistema

    n = size(U, 2);
    z = z(:);
    sol = zeros(n, 1);

    sol(n) = z(n) / U(n, n);

    for i = (n - 1):-1:1
        sol(i) = (z(i) - U(i, i + 1:n) * sol(i + 1:n)) / U(i, i);
    end
end