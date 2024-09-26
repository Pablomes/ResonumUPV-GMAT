function E2 = ErrorDiscreto(p, xi, fi)
% E2 = ErrorDiscreto(p, xi, fi)
% Calcula el error cuadratico de una aproximacion de puntos discretos
% PARAMETROS:
% p -> aproximacion de los puntos
% xi -> lista de coordenadas x
% fi -> lista de coordenadas y

xi = xi(:); fi = fi(:); n = length(fi);
E2 = 0; syms x;

for i = 1:n
    E2 = E2 + (subs(p, x, xi(i)) - fi(i)) * (subs(p, x, xi(i)) - fi(i));
end
end