function show_bar_custom(error_list, title_label)
    % SHOW_BAR_CUSTOM
    % Grafica un gráfico de barras agrupadas para cada lote
    % error_list: celda { {datos, label, style, color, lineSize, lineWidth}, ... }

    figure('Units', 'centimeters', 'Position', [2, 2, 14, 10]);
    hold on;

    legends = cell(1, numel(error_list));
    num_groups = numel(error_list{1}{1}); % Número de lotes
    num_series = numel(error_list);

    % Calcular posiciones para barras agrupadas
    bar_width = 0.8 / num_series; % ancho de cada barra
    x = 1:num_groups;

    for i = 1:num_series
        R = error_list{i}{1}; % datos (vector de errores por lote)
        label = error_list{i}{2};
        color = error_list{i}{4};

        % desplazamiento horizontal para agrupar
        offset = (i - (num_series + 1) / 2) * bar_width;

        % graficar barras
        bar(x + offset, R, bar_width, 'FaceColor', color, 'EdgeColor', 'none');

        legends{i} = label;
    end

    % Estética
    title(title_label, 'Interpreter', 'latex', 'FontSize', 12);
    xlabel('Snapshot $j$', 'Interpreter', 'latex', 'FontSize', 10);
    ylabel('Error', ...
           'Interpreter', 'latex', 'FontSize', 10);

    % Ajustar ticks solo a los lotes existentes
    xticks(x);
    xticklabels(arrayfun(@(j) sprintf('%d', j), 1:num_groups, 'UniformOutput', false));

    legend(legends, 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'northeast');
    set(gca, 'YScale', 'log'); % escala logarítmica
    set(gca, 'FontSize', 10, 'LineWidth', 1);

    % Guardar .eps
    clean_title = regexprep(title_label, '[^a-zA-Z0-9]', '');
    print('-depsc2', '-r300', fullfile('test', [clean_title, '_bars.eps']));

    hold off;
end
