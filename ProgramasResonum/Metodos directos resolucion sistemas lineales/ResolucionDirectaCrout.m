function sol = ResolucionDirectaCrout(diag, lDiag, uDiag, b)
% sol = ResolucionDirectaCrout(diag, lDiag, uDiag, b)
% Utiliza metodo de Crout para resolver el sistema Ax = b con
% matriz A tridiagonal
% PARAMETROS:
% diag -> diagonal principal de la matriz del sistema
% lDiag -> diagonal inferior a la principal
% uDiag -> diagonal superior a la principal
% b -> Vector independiente del sistema

    n = length(diag);
    diag = diag(:); lDiag = lDiag(:); uDiag = uDiag(:); b = b(:);

    % OBTENCION DE L Y U

    d = zeros(n, 1); u = zeros(n - 1, 1); l = zeros(n - 1, 1);

    d(1) = diag(1);
    u(1) = uDiag(1) / d(1);

    for i = 2:n
        l(i - 1) = lDiag(i - 1);
        d(i) = diag(i) - l(i - 1)*u(i - 1);
        
        if i < n
            u(i) = uDiag(i) / d(i);
        end
    end

    % SOLUCION DE Lz = d POR SUSTITUCION DIRECTA

    z = zeros(n, 1);

    z(1) = b(1) / d(1);
    
    for i = 2:n
        z(i) = (b(i) - l(i - 1)*z(i - 1)) / d(i);
    end

    % SOLUCION DE Ux = z POR SUSTITUCION INVERSA

    x = zeros(n, 1);

    x(n) = z(n);

    for i = (n - 1):-1:1
        x(i) = z(i) - u(i) * x(i + 1);
    end

    sol = x(:);

end