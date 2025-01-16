function I = Prob6MC(n)

    F = @(x, y, z) 2 * z * exp(-x^2);
    a = 0; b = 1;
    e = 1; f = 4; 

    sum = 0;
    i = 1;

    while i <= n
        xn = rand();
        yn = rand();
        zn = e + (f - e) * rand();

        if (yn <= xn)
            sum = sum + F(xn, yn, zn);
            i = i + 1;
        end
    end

    I = 1.5 * (1 / n) * sum;
end