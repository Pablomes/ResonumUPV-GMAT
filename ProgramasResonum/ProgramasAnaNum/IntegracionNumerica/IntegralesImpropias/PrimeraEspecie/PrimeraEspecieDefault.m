function I = PrimeraEspecieDefault(f, a, b)

I = quadgk(matlabFunction(f), a, b);

end