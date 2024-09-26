function x = Cramer(A, b)
d = det(A);
s = size(A);
s = s(2);
x = sym(zeros(1, s));
for i = 1:s
    x(i) = det([A(: ,1 : i - 1), b, A(:, i + 1 : s)]) / d;
end
end