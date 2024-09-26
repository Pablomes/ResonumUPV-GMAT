function [yk, xk, T] = Problema2b(F, a, b, alpha, beta, N, t0, t1, tol, maxiter)
% [yk, xk] = DisparoNoLinealSecante(F, a, b, alpha, beta, N, t0, t1, tol, maxiter)
% Obtiene la solucion del PC no lineal usando el metodo de la Secante para EDOS
% con condiciones de Dirichlet
% PARAMETROS:
% F -> sistema de EDOS de primer orden. Debe ser una funcion anonima de la
%       forma F = @(x, Y) [... ; ...]
% [a, b] -> Dominio del PC
% alpha -> valor de y(a)
% beta -> valor de y(b)
% N -> numero de subintervalos
% t0 -> Primera estimacion inicial de t
% t1 -> Segunda estimacion inicial de t
% tol -> tolerancia para S(t)
% maxiter -> maximo numero de iteraciones

h = (b-a) / N; xk = a:h:b; iter = 0; S = tol + 1; T = [t0, t1];

[Yk, ~] = MetodoRK4Sistemas(F, a, b, [t0, alpha], N);      y0 = Yk(1, end); dy0 = Yk(2, end);
[Yk, ~] = MetodoRK4Sistemas(F, a, b, [t1, alpha], N);      y1 = Yk(1, end); dy1 = Yk(2, end);

while and(S > tol, iter < maxiter)
    t2 = t1 - (y1 - dy1 - beta) * (t1 - t0) / (y1 - dy1 - y0 + dy0);
    T = [T, t2];
    [Yk, ~] = MetodoRK4Sistemas(F, a, b, [t2, alpha], N);
    y0 = y1; dy0 = dy1; y1 = Yk(1, end); dy1 = Yk(2, end); S = abs(y1 - dy1 - beta);
    t0 = t1; t1 = t2;
end

if S > tol
    yk = []; disp("No ha convergido")
else
    yk = Yk(1, :);
    hold off
    plot(xk, yk, "o--")
    axis padded
    hold on
end
end