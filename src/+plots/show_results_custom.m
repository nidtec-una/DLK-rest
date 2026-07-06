function show_results_custom(R_list, title_label, tol, save_folder)
    % Configuración de figura para una sola columna (SIAM suele usar ~13-16cm de ancho)
    fig = figure('Units', 'centimeters', 'Position', [2, 2, 13, 10], 'Color', 'w');

    hold on;
    num_curves = numel(R_list);
    legends = cell(1, num_curves);

    % Definición de estilos de línea y marcadores para máxima distinción en B/N

    for i = 1:num_curves
        R = R_list{i}{1};
        label = R_list{i}{2};

        % Cálculo de color en escala de grises (de negro a gris medio)
        grayVal = 0.6 - (i - 1) * (0.6 / max(1, num_curves - 1));
        grayColor = [grayVal, grayVal, grayVal];

        % Estilo cíclico de línea y marcador
        line_style = R_list{i}{3}(1:end - 1);
        marker = R_list{i}{3}(end);

        semilogy(R, line_style, ...
                 'Color', grayColor, ...
                 'Marker', marker, ...
                 'MarkerFaceColor', 'none', ...
                 'MarkerEdgeColor', grayColor, ...
                 'MarkerSize', 5, ...
                 'LineWidth', 1.2, ...
                 'MarkerIndices', round(linspace(1, length(R), 10))); % Solo 10 marcadores para limpieza

        legends{i} = label;
    end

    % Título y etiquetas con tipografía estándar de LaTeX
    % title(['\textbf{', title_label, '}'], 'Interpreter', 'latex', 'FontSize', 11);
    xlabel('Number of cycles ($k$)', 'Interpreter', 'latex', 'FontSize', 10);
    ylabel('Residual norm $\|r_k\|$', 'Interpreter', 'latex', 'FontSize', 10);

    % Leyenda sin recuadro (típico de SIAM) o con recuadro fino
    lgd = legend(legends, 'Interpreter', 'latex', 'FontSize', 9, ...
                 'Location', 'northeast', 'EdgeColor', 'none', 'Color', 'none');

    % Ajustes de ejes (Ticks hacia afuera y fuente Times)
    ax = gca;
    set(ax, 'YScale', 'log', 'TickLabelInterpreter', 'latex', ...
        'FontSize', 10, 'LineWidth', 0.8, 'FontName', 'Times', ...
        'XMinorTick', 'on', 'YMinorTick', 'on');

    % Ajuste de límites
    ylim([tol, 1.2]); % Un pequeño margen superior sobre 1
    grid on;
    set(ax, 'GridLineStyle', ':', 'GridAlpha', 0.3);

    % Exportar en formato EPS (Vectorial para LaTeX)
    clean_title = regexprep(title_label, '[^a-zA-Z0-9]', '');
    % Exportar con resolución para publicación
    if ~exist(save_folder, 'dir')
        mkdir(save_folder);
    end
    exportgraphics(fig, fullfile(save_folder, [clean_title, '.eps']), 'ContentType', 'vector');

    hold off;
end
