function [RD, RE] = Rectangulos(f, a, b, n)
% [RD, RE] = Rectangulos(f, a, b, n)
% Aproxima la integral a base de rectangulos tanto en defecto como en
% exceso
% PARAMETROS:
% f -> funcion anonima a integrar
% [a, b] -> dominio
% n -> numero de rectangulos

RD = 0; RE = 0;
h = (b - a) / n;

for i = 0:(n-1)
    RD = RD + f(a + h * i);
    RE = RE + f(a + (h * (i + 1)));
end

RD = RD * h;
RE = RE * h;

end