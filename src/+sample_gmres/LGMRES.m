% --- Funcion de muestreo del LGMRES con parada automática sin flag ---
function [MatR, xf, vecNormRrel, cycles] = LGMRES(A, b, x0, m, l, tol, number_of_samples)

    % Inicialización
    R = b - A * x0;
    normR0 = norm(R);

    if normR0 == 0
        MatR = [];
        xf = x0;
        vecNormRrel = 0;
        cycles = 0;
        return
    end

    % Verificar number_of_samples > = m + l
    if number_of_samples < l + 1
        error('El numero de ciclos de muestreo debe ser mayor o igual al del #{subespacion de enrriquecimiento}+1 = l+1');
    end

    number_of_samples = number_of_samples + 1; % Ajuste por el residuo inicial
    [xf, vecNormRrel_t, MatR] = utils.lgmres_v2(A, b, m, l, tol, number_of_samples, x0);
    vecNormRrel = vecNormRrel_t';
    cycles = size(MatR, 2) - 1; % Numero de ciclos realizados
end
