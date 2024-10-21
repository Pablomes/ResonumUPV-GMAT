function [ACOC] = fACOC(incr)
% ACOC = fACOC(incr)
% Calcula el ACOC
% PARAMETROS:
% incr -> lista de incrementos del algoritmo

    ACOC = [];
    for i = 2 : (length(incr) - 1)
        ACOC = [ACOC, (log(incr(i + 1) / incr(i)) / log(incr(i) / incr(i - 1)))];
    end
end