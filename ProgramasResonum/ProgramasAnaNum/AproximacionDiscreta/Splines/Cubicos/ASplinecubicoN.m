function [S,p]= ASplinecubicoN(x,y) 
% [S,p]= SplinecubicoN(x,y)
% Spline cúbico Natural
% introduces los dos vectores de datos y devuelve la matriz de polinomios
% que tendrá cuatro columnas y un conjunto de datos que se puede 
% operar con ppval(p,valor)

%% Cálculo de la matriz A

% esto calcula las diferencias entre cada elemento del vector 
% si el primer elemento es 5 y el segundo 7, devuelve 2
% se usa para calcular longitudes de los segmentos entre puntos de control
% h es el vector que guarda los (x-xi)

h=diff(x); 

% Extracción de los coeficientes de c para obtener sus valores.
% se crea una matriz que en cada fila solo tiene tres elementos
% y van moviéndose horizontalmente
% recordemos que en el spline natural la primera fila son todo
% ceros y un uno, al principio y al final respectivamente

% esto lo escribimos así para usar la función Crout

% diagonal superior

ds=[0 h(2:end)];

% diagonal inferior

di=[h(1:end-1) 0];

% diagonales principales.

dp= [1 2*(h(1:end-1)+h(2:end)) 1];


%% Cálculo de los términos independientes

% aquí hacemos lo mismo que con x
% calculamos las diferencias. a(i+1)-a(i)

a=y;

l=diff(a);

% a partir de aquí vamos a crear el término independiente del sistema

% calcula las pendientes de las rectas usando la fórmula de diferencia
% dividida

m=l./h;

% Calcula las diferencias entre las pendientes calculadas consecutivas. 
% Estas diferencias son necesarias para calcular los términos 
% independientes del sistema de ecuaciones

M=diff(m);

% Calcula el vector de términos independientes del sistema de ecuaciones

v=3*[0 M(1:end) 0];

% Resuelve el sistema de ecuaciones lineales utilizando el método de 
% Crout para obtener los coeficientes del spline cúbico

c=Crout(dp,ds,di,v);

% Calcula los coeficientes d para cada segmento del spline cúbico. 
% Estos coeficientes se calculan a partir de los coeficientes c 
% obtenidos anteriormente y las longitudes de los intervalos entre 
% puntos de control h. 
% Fórmula en los apuntes

d=[];
for i=1:length(x)-1
    di=(c(i+1)-c(i))/(3*h(i));
    d=[d di];
end

% Lo mismo que con d

b=[];
for i=1:length(x)-1
    bi=(a(i+1)-a(i))/h(i)-h(i)/3*(2*c(i)+c(i+1));
    b=[b bi];
end

%% Ajustar las dimensiones de los vectores a, b, c y d y creación de S y p

c=c(1:end-1);
a=a(1:end-1);
a=a(:);
b=b(:);
c=c(:);
d=d(:);

% importante que el orden sea [d c b a]

S=[d c b a]; 

% Crea la pieza de polinomios utilizando los puntos de control x y los 
% coeficientes S calculados.

p=mkpp(x,S);

end

% bendito ChatGPT