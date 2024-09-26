function n = Problema8_Tema4(f, a, b, x0, tol)
    e = 1; n = -1; syms x;

    while e > tol
        n = n + 1;
        p = AproximacionPade(f, n, n, a, b);
        e = abs(subs(f, x, x0) - subs(p, x, x0));
    end
end