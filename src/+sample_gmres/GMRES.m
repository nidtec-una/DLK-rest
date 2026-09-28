% --- Funcion de muestreo del GMRES con parada automática ---
% OBS: El almacenamiendo de r0/|r0| no cuenta como un ciclo
function [MatR, xf, vecNormRrel, cycles] = GMRES(A, b, x0, m, tol, number_of_samples, extra_info)
    if nargin < 7
        extra_info = [];
    end

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

        [~, xf] = evalc('gmres(A, b, m, tol, 1, [], [], x0)');
        R = b - A * xf;
        nR = norm(R);
        rRel = nR / max(normb, eps);

        % Guarda resultados
        vecNormRrel(i) = rRel;
        MatR(:, i) = R / max(nR, eps);

        % Progress loading indicator
        curr_c = i - 1;
        if vecNormRrel(end) < tol || i > number_of_samples
            pct = 100;
        else
            pct = min(99, round((curr_c / number_of_samples) * 100));
        end
        nb = round((pct / 100) * 20);
        bar_str = [repmat('=', 1, nb), repmat('-', 1, 20 - nb)];

        if isstruct(extra_info)
            curr_tot_cycles = extra_info.tot_cycles + curr_c;
            fprintf(1, '\r  GMRES Progress: [%s] %3d%% (cycle %d/%d) (Snapshot %d, m = %d | cycles %d/%d)', ...
                    bar_str, pct, curr_c, number_of_samples, extra_info.j, extra_info.m, curr_tot_cycles, extra_info.cycleMax);
        else
            fprintf(1, '\r  GMRES Progress: [%s] %3d%% (cycle %d/%d)', bar_str, pct, curr_c, number_of_samples);
        end
        drawnow;

        x0 = xf;
        i = i + 1;
    end
    fprintf(1, '\n');
    cycles = size(MatR, 2) - 1;
    xf = x0;
end
