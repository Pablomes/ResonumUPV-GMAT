function I = Prob1ExamenPasado(F, n)
    
    [xk, yk] = Puntos(n);

    sum = 0;

    for i = 1:n
        sum = sum + F(xk(i), yk(i));
    end

    I = 12 * (1 / n) * sum;
end

function [xn, yn] = Puntos(n)
    
    xn = zeros(n, 1);
    yn = zeros(n, 1);

    i = 1;

    while i <= n
        xk = (rand() - 0.5) * 4;
        yk = (rand() - 0.5) * 4;

        if (xk <= -1 || xk >= 1 || yk <= -1 || yk >= 1)
            xn(i) = xk;
            yn(i) = yk;

            i = i + 1;
        end
    end
end