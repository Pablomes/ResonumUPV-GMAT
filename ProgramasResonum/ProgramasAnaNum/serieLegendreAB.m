function SL = serieLegendreAB(f,N,a,b)
% SL = serieLegendreAB(f,N,a,b) obtiene los N primeros 
% términos de la serie de Legendre de f(x) en un 
% intervalo [a,b], siendo f una función simbólica
syms x t
% X=(b-a)/2*t+(b+a)/2
ft=subs(f,x,(b-a)/2*t+(b+a)/2);
SL=0;

for k=0:N
    P = legendreP(k,t);
    C = (2*k+1)/2*int(ft*P,t,-1,1); % primero lo hacemos entre -1 y 1
    C = double(C);
    SL = SL+ C*P;
end
% ahora deshacemos el cambio de variable
% t=(2*x-b-a)/(b-a)
SL = subs(SL,t,(2*x-b-a)/(b-a));