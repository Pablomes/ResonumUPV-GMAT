function Tasa(incr)
% Tasa(incr)
% Calcula y muestra las tasas de error
% PARAMETROS:
% incr -> lista de incrementos del algoritmo

    k = 1:length(incr) - 2;
    figure,
    plot(k, incr(2:end - 1)./incr(1:end - 2), "o--")
    hold on
    plot(k, incr(2:end - 1)./incr(1:end - 2).^2, "o--")
    plot(k, incr(2:end - 1)./incr(1:end - 2).^3, "o--")
    legend("Lineal", "Cuadratica", "Cubica")
    xlabel("k"),grid
    axis([k(1), k(end), 0, 4])
end