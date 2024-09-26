function SC = serieChebyshevAB(f,N,a,b)
% SC = serieChebyshev(f,N,a,b) obtiene los N primeros 
% términos de la serie de Chebyshev de f(x) en un 
% intervalo [a,b], siendo f una función simbólica
syms x t
w=1/sqrt(1-t^2);
ft=subs(f,x,(b-a)/2*t+(b+a)/2);
P=chebyshevT(0,t);
a0=2/pi*int(ft*P*w,t,-1,1);
SC=a0/2;

for k=1:N
    P=chebyshevT(0,t);
    ak=2/pi*int(ft*P*w,t,-1,1);
    SC=SC+ak*P;
end
SC= subs(SC,t,(2*x-b-a)/(b-a));
