function [yk, xk] = Problema3(N, t0, t1, tol, maxiter)
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

F = @(x, Y) [Y(2); Y(3); -Y(3)]; a = 1; b = 2; alpha1 = 0; alpha3 = 1; beta = 2.69315;

h = (b-a) / N; xk = a:h:b; iter = 0; S = tol + 1;

[Yk, ~] = MetodoRK4Sistemas(F, a, b, [alpha1; t0; alpha3], N); dya0 = Yk(2, 1); dyb0 = Yk(2, end);
[Yk, ~] = MetodoRK4Sistemas(F, a, b, [alpha1; t1; alpha3], N); dya1 = Yk(2, 1); dyb1 = Yk(2, end);

while and(S > tol, iter < maxiter)
    t2 = t1 - (dya1 + dyb1 - beta) * (t1 - t0) / (dya1 + dyb1 - dya0 - dyb0);
    [Yk, ~] = MetodoRK4Sistemas(F, a, b, [alpha1; t2; alpha3], N);
    dya0 = dya1; dyb0 = dyb1; dya1 = Yk(2, 1); dyb1 = Yk(2, end);
    S = abs(dya1 + dyb1 - beta);
    t0 = t1; t1 = t2; iter = iter + 1;
end

if S > tol
    yk = []; disp("No ha convergido")
else
    %yk = Yk(1, :);
    yk = Yk;
    hold off
    plot(xk, Yk(1, :), "o--")
    axis padded
    hold on
end
end