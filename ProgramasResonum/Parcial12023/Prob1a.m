function [A, b] = Prob1a(n)

    A = diag((1:n) * 10);
    b = zeros(n, 1);

    for i = 1:n
        for j = (i + 1):n
            A(i, j) = i - j + 1;
        end

        b(i) = (-1)^(i + 1);
    end
end