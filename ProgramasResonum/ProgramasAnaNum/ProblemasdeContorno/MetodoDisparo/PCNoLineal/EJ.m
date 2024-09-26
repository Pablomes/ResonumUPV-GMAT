function [yk, xk] = EJ(F, a, b, N, t0, t1, tol, maxiter)

h = (b-a) / N; xk = a:h:b; iter = 0; S = tol + 1;

[Yk, ~] = MetodoRK4Sistemas(F, a, b, [17; t0; 34], N);  y0a = Yk(2, 1);    y0b = Yk(2, end);
[Yk, ~] = MetodoRK4Sistemas(F, a, b, [17; t1; 34], N);  y1a = Yk(2, 1);    y1b = Yk(2, end);

while and(S > tol, iter < maxiter)
    t2 = t1 - (y1a + y1b + 7) * (t1 - t0) / (y1a + y1b - y0a - y0b);
    [Yk, ~] = MetodoRK4Sistemas(F, a, b, [17; t2; 34], N);
    y0a = y1a; y0b = y1b; y1a = Yk(2, 1); y1b = Yk(2, end); S = abs(y1a + y1b + 7);
    t0 = t1; t1 = t2; iter = iter + 1;
end

if S > tol
    yk = []; disp("No ha convergido")
else
    yk = Yk(1, :);
    hold off
    plot(xk, yk, "o--")
    axis padded
    hold on
end
end