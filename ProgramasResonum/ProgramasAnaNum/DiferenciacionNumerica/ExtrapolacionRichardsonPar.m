function dVal = ExtrapolacionRichardsonPar(N, H, decExact, x0)
% dVal = ExtrapolacionRichardsonPar(N, H, decExact, x0)
% Aproxima el diferencial a base de la extrapolacion de Richardson para las
% potencias pares de h
% PARAMETROS:
% N -> Aproximacion N(h). Debe darse como una expresion simbolica que
%       contenga h, y x para representar x0.
% H -> Valor de la separacion
% decExact -> numero de decimales exactos de la solucion
% x0 -> punto en el que se quiere la derivada

syms x h; k = -(decExact) / (log(H) / log(10));
if ceil(k) == k
    if mod(ceil(k), 2) == 0
        k = ceil(k) + 2;
    else
        k = ceil(k) + 1;
    end
else
    k = ceil(k);
    if mod(k, 2) == 1
        k = k + 1;
    end
end

F = sym(zeros(k, k));N = subs(N, x, x0);
for row = 1:k
    F(row, 1) = subs(N, h, H / (2^(row - 1)));
end

for col = 2:k
    for row = 1:(k + 1 - col)
        F(row, col) = (1 / (4 ^(col - 1) - 1)) * ((4 ^ (col - 1)) * F(row + 1, col - 1) - F(row, col - 1));
    end
end

dVal = vpa(F(1, k));

end