% --- Funcion de muestreo del GMRES con parada automática ---
% OBS: El almacenamiendo de r0/|r0| no cuenta como un ciclo
function [MatR, xf, vecNormRrel, cycles] = GMRES(A, b, x0, m, tol, number_of_samples)

    % Inicialización
    MatR = [];
    vecNormRrel = [];
    R = b - A * x0;
    normR0 = norm(R);

    if normR0 == 0
        MatR = [];
        xf = x0;
        vecNormRrel = 0;
        cycles = 0;
        return
    end

    normb = norm(b);
    vecNormRrel(1) = normR0 / max(normb, eps);
    MatR(:, 1) = R / max(normR0, eps);

    i = 2;
    while i <= number_of_samples + 1 && vecNormRrel(end) >= tol

        xf = gmres(A, b, m, tol, 1, [], [], x0);
        R = b - A * xf;
        nR = norm(R);
        rRel = nR / max(normb, eps);

        % Guarda resultados
        vecNormRrel(i) = rRel;
        MatR(:, i) = R / max(nR, eps);

        x0 = xf;
        i = i + 1;
    end
    cycles = size(MatR, 2) - 1;
    xf = x0;
end
