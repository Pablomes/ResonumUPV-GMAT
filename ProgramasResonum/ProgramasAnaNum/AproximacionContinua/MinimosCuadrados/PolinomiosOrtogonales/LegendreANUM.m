function [ Pn ] = LegendreANUM (n , var )
% [Pn] = legendreANUM (n, var )
% obtiene el polinomio de Legendre Pn( var ) de grado n
% PARAMETROS:
% n -> grado
% var -> expresion simbolica (suele ser x)

pn {1}=1;
if n >0
pn {2}= var ;
if n >1
for k =3: n +1
bk = int ( var .* pn {k -1}^2 , var , -1 ,1) /...
int ( pn {k -1}.^2 , var , -1 ,1) ;
ck = int ( var .* pn {k -1}.* pn {k -2} , var , -1 ,1) /...
int ( pn {k -2}^2 , var , -1 ,1) ;
pn {k }= simplify (( var - bk )* pn {k -1} - ck * pn {k -2}) ;
end
end
end
pn = pn (:) ; Pn = pn { end };
