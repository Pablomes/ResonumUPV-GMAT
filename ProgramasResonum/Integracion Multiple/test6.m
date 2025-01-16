function I=test6(f,a,b,c, d, e,g,n)
%Definir area
Volumen=1.5;
k=1;
x=zeros(n,1);
y=zeros(n,1);
z=zeros(n,1);
while k<n
    u=rand(1,1);
    v=rand(1,1);
    w=rand(1,1);
     %Cambiar dominio con x=a+(b-a)*t
    u=a+(b-a)*u;
    v=c+(d-c)*v;
    w=e +(g-e)*w;
    %Insertar condicion
    cond1=1;
    cond2=v<u;
    cond3=1;
    if and(cond1,and(cond2,cond3))
        %plot(u,v,"*");
        %hold on
        x(k)=u;
        y(k)=v;
        z(k)=w;
        k=k+1;
        
    end
    
end

S=feval(f,x,y,z);
I=(1/n)*Volumen*sum(S);
plot3(x,y,z,"*");