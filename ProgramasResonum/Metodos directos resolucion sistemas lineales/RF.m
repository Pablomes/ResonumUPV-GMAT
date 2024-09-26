function RF(A, b)
% RF(A, b)
% Usa Rouche-Frobenius para identificar el sistema y lo muestra por pantalla
% PARAMETROS:
% A -> Matriz cuadrada del sistema
% b -> Vector del sistema

    n = size(A, 2);
    b = b(:);

    Am = [A, b];

    r = rank(A);
    rAm = rank(Am);

    if (r == rAm)
        if (r == n)
            disp("Es sistema compatible determinado.")
        else 
            disp("Es sistema compatible indeterminado.")
        end
    else
        disp("Es sistema indeterminado.")
    end

end