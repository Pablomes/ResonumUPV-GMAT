function [N, p, E] = Problema2(f, n, tol)
N = 0; E = 1; syms x;

while E > tol
    hold off
    p = AproximacionLegendrePGorro(f, N);
    
    xi = linspace(-1, 1, n + 1);

    for i = 1:length(xi)
        e = abs(subs(f, x, xi(i)) - subs(p, x, xi(i)));
        if or(e>E, i == 1)
            E = e;
        end
    end

    if E > tol
        N = N +1;
    end
end
end