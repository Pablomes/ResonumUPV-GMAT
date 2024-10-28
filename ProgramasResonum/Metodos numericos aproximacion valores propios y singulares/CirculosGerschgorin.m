function [D] = CirculosGerschgorin(A)
% [D] = CirculosGerschgorin(A)
% Calcula y devuelve una tabla con los circulos de Gerschgorin en el
% formato D(i) = [Centro, Radio]
% PARAMETROS:
% A -> Matriz A de la que se quieren los circulos

    n = size(A, 1);
    D = zeros(n, 2);

    for i = 1:n
        D(i, 1) = A(i, i);
        sum = 0;
        for j = 1:n
            if i~=j
                sum = sum + abs(A(i, j));
            end
        end

        D(i, 2) = sum;
    end
end