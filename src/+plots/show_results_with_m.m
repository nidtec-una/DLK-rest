function show_results_with_m(R_list, title_label, m_limits, m_values, p_values, save_folder)
    % Versión SIAM con desplazamiento de gradiente y alineación exacta

    % -----------------------------
    % Configurar p_values si no viene
    % -----------------------------
    if nargin < 5 || isempty(p_values)
        p_values = 5 * ones(size(m_values));
    end

    start_cycle = [0, cumsum(p_values(1:end - 1))];
    end_cycle_snapshots = cumsum(p_values);
    maxX = length(R_list{1}{1});
    % xEnd = 36;
    xEnd = 43;

    end_cycle_snapshots;

    % Proporciones de figura optimizadas
    fig = figure('Units', 'centimeters', 'Position', [2, 2, 17, 18], 'Color', 'w');

    darkGray = [0.15, 0.15, 0.15];
    midGray  = [0.45, 0.45, 0.45];

    xticks_vec = 0:5:xEnd;

    %% ============================================================
    % (1) CURVAS DE RESIDUOS
    % ============================================================
    subplot(3, 1, 1);
    hold on;

    num_series = numel(R_list);

    for i = 1:num_series
        R = R_list{i}{1};
        label = R_list{i}{3};

        % Color dinámico (Escala de grises SIAM)
        % gVal = (num-1) * (0.5 / max(1, 4-1))
        % c = [gVal, gVal, gVal];
        c = [0.6 0.6 0.6] - R_list{i}{5};

        if numel(R) > maxX + 1
            R = R(1:maxX + 1);
        end

        line_style = R_list{i}{4}(1:end - 1);
        marker = R_list{i}{4}(end);

        semilogy(0:(numel(R) - 1), R, line_style, ...
                 'Color', c, ...
                 'Marker', marker, ...
                 'MarkerFaceColor', 'none', ...
                 'MarkerEdgeColor', c, ...
                 'MarkerSize', 4, ...
                 'LineWidth', 1, ...
                 'MarkerIndices', round(linspace(1, length(R), 12)));
    end

    % title(['\textbf{', title_label, '}'], 'Interpreter', 'latex', 'FontSize', 11);
    ylim([1e-13, max(R_list{1}{1})]);
    ylabel('$\|r_k\|/\|r_0\|$', 'Interpreter', 'latex', 'FontSize', 11);
    xlim([0, xEnd]);
    ylim([1e-10, 1]);
    ax = findobj(gcf, 'Type', 'axes');
    set(ax, 'XTick', xticks_vec);
    set(gca, 'YScale', 'log', 'FontSize', 10, 'LineWidth', 0.8, 'TickLabelInterpreter', 'latex');
    grid on;

    %% ============================================================
    % (2) BARRAS DE ERROR (GRADIENTE DESPLAZADO)
    % ============================================================
    subplot(3, 1, 2);
    hold on;

    % Recolección de errores para escala
    all_errors = [];
    for i = 1:num_series
        all_errors = [all_errors, R_list{i}{2}(:)'];
    end
    min_err = min(all_errors);
    max_err = max(all_errors);

    % Mapa de colores (Blanco a Gris Oscuro según tu estilo original)
    n_colors = 256;
    cmap = repmat(linspace(0.8, 0.2, n_colors)', 1, 3);

    % --- INTEGRACIÓN DEL DESPLAZAMIENTO Y EXTENSIÓN FINAL ---
    gradient_shift = 0.02;

    for i = 1:num_series
        R_error = R_list{i}{2};
        num_snap = min(numel(R_error), numel(start_cycle));

        for j = 1:num_snap
            x1 = start_cycle(j);

            % Si es la última barra del ciclo, extender hasta xEnd
            % if j == num_snap
            %    if end_cycle_snapshots(j) < maxX
            %        x2 = maxX-1;
            %    else
            %        x2 = end_cycle_snapshots(j);
            %    end
            % else
            %    x2 = end_cycle_snapshots(j);
            % end
            x2 = end_cycle_snapshots(j);

            y = R_error(j);

            idx = round((y - min_err) / (max_err - min_err) * (n_colors - 1)) + 1;
            idx = max(1, min(n_colors, idx));

            % Dibujamos la barra con el shift incorporado
            fill([x1 x2 x2 x1] + gradient_shift, [min_err * 0.9, min_err * 0.9, y, y], ...
                 cmap(idx, :), 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 0.2, 'FaceAlpha', 0.9);
        end
    end

    colormap(gca, cmap);
    caxis([min_err, max_err]);

    % --- COLORBAR CON DESPLAZAMIENTO A LA DERECHA ---
    cb = colorbar('eastoutside');
    cb.TickLabelInterpreter = 'latex';
    % cb.Label.String = 'Relative Error $e$';
    % cb.Label.Interpreter = 'latex';

    cb_pos = cb.Position;
    cb_pos(1) = cb_pos(1) + 0.08; % Mover colorbar a la derecha para no amontonar
    cb.Position = cb_pos;

    ylabel('$e = \|\hat{S} - S\|_F/\|S\|_F$', 'Interpreter', 'latex', 'FontSize', 11);
    xlim([0, xEnd]);
    ax = findobj(gcf, 'Type', 'axes');
    set(ax, 'XTick', xticks_vec);
    set(gca, 'YScale', 'log', 'FontSize', 10, 'LineWidth', 0.8, 'TickLabelInterpreter', 'latex');
    grid on;

    %% ============================================================
    % (3) VALORES DE m
    % ============================================================
    subplot(3, 1, 3);

    m_x = repelem(m_values(:).', p_values(:).');
    if numel(m_x) < maxX + 1
        m_x = [m_x, repmat(m_x(end), 1, maxX + 1 - numel(m_x))];
    else
        m_x = m_x(1:maxX + 1);
    end

    stairs(0:maxX, m_x, 'Color', darkGray, 'LineWidth', 1.2);

    % yline(m_limits(1), ':', 'Color', midGray, 'Label', ['$m_{min} = ', num2str(m_limits(1)), '$'], ...
    %    'Interpreter', 'latex', 'LabelHorizontalAlignment','left');
    yline(m_limits(2), '--', 'Color', midGray, 'Label', ['$m_{max} = ', num2str(m_limits(2)), '$'], ...
          'Interpreter', 'latex', 'LabelHorizontalAlignment', 'left');

    xlabel('Number of cycles ($k$)', 'Interpreter', 'latex', 'FontSize', 11);
    ylabel('$m$', 'Interpreter', 'latex', 'FontSize', 11);
    xlim([0, xEnd]);
    ylim([m_limits(1) * 0.8, m_limits(2) * 1.2]);
    ax = findobj(gcf, 'Type', 'axes');
    set(ax, 'XTick', xticks_vec);
    set(gca, 'FontSize', 11, 'LineWidth', 0.8, 'TickLabelInterpreter', 'latex');
    grid on;

    % Sincronización de ejes
    linkaxes(findobj(gcf, 'Type', 'axes'), 'x');

    %% ============================================================
    % GUARDAR EPS (VECTORIAL)
    % ============================================================
    clean_title = regexprep(title_label, '[^a-zA-Z0-9]', '');
    if ~exist(save_folder, 'dir')
        mkdir(save_folder);
    end
    exportgraphics(fig, fullfile(save_folder, [clean_title, '_with_m.eps']), 'ContentType', 'vector');

end
