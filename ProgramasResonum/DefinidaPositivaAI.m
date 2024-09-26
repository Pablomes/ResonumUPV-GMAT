function pos = DefinidaPositivaAI(A)
% pos = DefinidaPositiva(A)
% Determina si la matriz es definida positiva, devuelve 1 si lo es, 0 si no
% PARAMETROS:
% A -> Matriz cuadrada

    d = A(1, 1);
    k = 2;

    while d > 0 && k <= size(A, 2)
        d = det(A(1:k, 1:k));
        k = k + 1;
    end

    if d > 0
        disp("DEFINIDA POSITIVA.")
        pos = 1;
    else
        disp("NO DEFINIDA POSITIVA.")
        pos = 0;
    end

    
end