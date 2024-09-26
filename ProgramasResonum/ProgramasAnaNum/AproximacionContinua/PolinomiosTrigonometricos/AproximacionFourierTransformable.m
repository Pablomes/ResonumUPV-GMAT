function p = AproximacionFourierTransformable(f, n, P, c)
% p = AproximacionFourierPeriodo(f, n, P)
% Aproxima la funcion usando la serie de Fourier de grado n
% PARAMETROS:
% f -> expresión simbolica a aproximar
% n -> grado
% P -> periodo de aproximacion
% c -> centro de aproximacion

syms x;
p = 0; w = 2 * pi / P;

a = (2 / P) * int(f, x, -P / 2 + c, P / 2 + c);
p = p + (a / 2);

for k = 1:(n)
    a = (2 / P) * int(f * cos(k * x * w), x, -P / 2 + c, P / 2 + c);
    b = (2 / P) * int(f * sin(k * x * w), x, -P / 2 + c, P / 2 + c);

    p = p + (a * cos(k * x * w) + b * sin(k * x * w));
end
%if n > 0
%    a = (1 / pi) * int(f * cos(n * x), x, -pi, pi);
%    p = p + (a * cos(n * x));
%end

fplot(f, [-P/2 + c, P/2 + c])
hold on
fplot(p, [-P/2 + c, P/2 + c])
legend("f(x)", "p" + n + "(x)")

end