function [yk, xk] = test(F,Fy1,Fy2,a,b,alpha,beta,N,t0,tol,maxiter)
% [yk, xk] = DisparoNoLinealNewton(F, Fy1, Fy2, a, b, alpha, beta, N, t0, tol, maxiter)
% Obtiene la solucion del PC no lineal usando el metodo de Newton para EDOS
% de la forma y''(x) = F(x, y, y') con condiciones de Dirichlet
% PARAMETROS:
% F -> funcion F(x, y, y') de la EDO. Debe ser una funcion anonima de la
%       forma F = @(x, y, dy) ...
% Fy1 -> derivada respecto de y de la funcion F(x, y, y'). Debe ser una 
%       funcion anonima de la forma Fy1 = @(x, y, dy) ...
% Fy2 -> derivada respecto de y' de la funcion F(x, y, y'). Debe ser una 
%       funcion anonima de la forma Fy2 = @(x, y, dy) ...
% [a, b] -> Dominio del PC
% alpha -> valor de y(a)
% beta -> valor de y(b)
% N -> numero de subintervalos
% t0 -> Estimacion inicial de t
% tol -> tolerancia para S(t)
% maxiter -> maximo numero de iteraciones

h = (b - a) / N; xk = a:h:b; iter = 0; S = tol + 1;
FN = @(x, Y) [Y(2); F(x, Y(1), Y(2)); Y(4); Fy1(x, Y(1), Y(2)) * Y(3) + Fy2(x, Y(1), Y(2)) * Y(4)];

[Yk, ~] = MetodoRK4Sistemas(FN, a, b, [t0;alpha;1;0], N); 
y1b = Yk(1, end); dy1b = Yk(2, end); y3b = Yk(3, end); dy3b = Yk(4, end);

while and(S > tol, iter < maxiter)

    t1 = t0 - (dy1b - beta) / dy3b;
    [Yk, ~] = MetodoRK4Sistemas(FN, a, b, [t1;alpha;1;0], N); 
    y1b = Yk(1, end); dy1b = Yk(2, end); y3b = Yk(3, end); dy3b = Yk(4, end); S = abs(dy1b - beta);
    t0 = t1; iter = iter + 1;
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

yk = Yk(1, :);

end
