function p = AproximacionFourier(f, n)
% p = AproximacionFourier(f, n)
% Aproxima la funcion usando la serie de Fourier de grado n
% PARAMETROS:
% f -> expresión simbolica a aproximar
% n -> grado

syms x;
p = 0;

a = (1 / pi) * int(f, x, -pi, pi);
p = p + (a / 2);

for k = 1:(n)
    a = (1 / pi) * int(f * cos(k * x), x, -pi, pi);
    b = (1 / pi) * int(f * sin(k * x), x, -pi, pi);

    p = p + (a * cos(k * x) + b * sin(k * x));
end
%if n > 0
%    a = (1 / pi) * int(f * cos(n * x), x, -pi, pi);
%    p = p + (a * cos(n * x));
%end

fplot(f, [-pi, pi])
hold on
fplot(p, [-pi, pi])
legend("f(x)", "p" + n + "(x)")

end