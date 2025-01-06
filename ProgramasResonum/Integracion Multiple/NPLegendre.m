function [Nodos,Pesos]=NPLegendre(n)
% [Nodos, Pesos] = NPLegendre(n)
% Calcula los nodos y los pesos del polinomio de Legendre dada una n
% PARAMETROS:
% n -> numero de nodos
    r=1;
    q=[1 0];
    for k=1:n-1
        p=(1/(k+1))*((2*k+1)*[q 0]-k*[0 0 r]);
        r=q;
        q=p;
    end

    Nodos=roots(p);

    dp=polyder(p);
    z=polyval(dp,Nodos);
    Pesos=2./((ones(n,1)-Nodos.^2).*(z.^2));
end