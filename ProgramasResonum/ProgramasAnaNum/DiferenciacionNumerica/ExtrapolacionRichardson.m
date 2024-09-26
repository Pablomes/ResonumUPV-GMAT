function dVal = ExtrapolacionRichardson(N, H, decExact, x0)
% dVal = ExtrapolacionRichardson(N, H, decExact, x0)
% Aproxima el diferencial a base de la extrapolacion de Richardson para
% todas las potencias de h
% PARAMETROS:
% N -> Aproximacion N(h). Debe darse como una expresion simbolica que
%       contenga h, y x para representar x0.
% H -> Valor de la separacion
% decExact -> numero de decimales exactos de la solucion
% x0 -> punto en el que se quiere la derivada

syms x h; k = -(decExact) / (log(H) / log(10));
if ceil(k) == k
    k = ceil(k) + 1;
else
    k = ceil(k);
end

F = sym(zeros(k, k));N = subs(N, x, x0);
for row = 1:k
    F(row, 1) = subs(N, h, H / (2^(row - 1)));
end

for col = 2:k
    for row = 1:(k + 1 - col)
        F(row, col) = (1 / (2 ^(col - 1) - 1)) * ((2 ^ (col - 1)) * F(row + 1, col - 1) - F(row, col - 1));
    end
end

dVal = vpa(F(1, k));

end