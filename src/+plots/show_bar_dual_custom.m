function show_bar_dual_custom(data_list, title_label, y_label, log_scale, save_folder, x_labels)

    if nargin < 5 || isempty(save_folder)
        save_folder = pwd;
    end
    if ~exist(save_folder, 'dir')
        mkdir(save_folder);
    end

    % --- 1. FILTRADO DE DATOS (Eliminar NC/NaN) ---
    num_series = numel(data_list);
    num_total_matrices = numel(data_list{1}{1});
    valid_mask = true(1, num_total_matrices);
    for i = 1:num_series
        valid_mask = valid_mask & ~isnan(data_list{i}{1})';
    end

    filtered_data = cell(1, num_series);
    for i = 1:num_series
        filtered_data{i} = data_list{i}{1}(valid_mask);
    end
    active_ids = (1:num_total_matrices);
    active_ids = active_ids(valid_mask);

    % --- 2. FIGURA - Proporciones SIAM ---
    fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 16 8.5]);
    tlo = tiledlayout(1, 1, 'Padding', 'loose', 'TileSpacing', 'compact');
    ax = nexttile(tlo);
    hold(ax, 'on');

    % --- 3. DIBUJO DE BARRAS ---
    num_active_groups = numel(active_ids);
    bar_width = 0.8 / num_series;
    x_pos = 1:num_active_groups;
    handles = gobjects(1, num_series);
    legends = cell(1, num_series);

    for i = 1:num_series
        R = filtered_data{i};
        offset = (i - (num_series + 1) / 2) * bar_width;

        handles(i) = bar(ax, x_pos + offset, R, bar_width, ...
                         'FaceColor', data_list{i}{3}, ...
                         'EdgeColor', [0.15 0.15 0.15], ...
                         'LineWidth', 0.5, ...
                         'FaceAlpha', 0.85);

        legends{i} = data_list{i}{2};
    end

    % --- 4. AJUSTE DE LIMITES (PADDING TOTAL) ---
    all_values = cell2mat(filtered_data');
    max_val = max(all_values(:));
    min_val_all = min(all_values(:));
    min_pos_val = min(all_values(all_values > 0));

    % PADDING VERTICAL (Y)
    if log_scale
        set(ax, 'YScale', 'log');
        top_lim = 10^(log10(max_val) + 1.5);
        bottom_lim = 10^(log10(min_pos_val) - 0.5);
        ylim(ax, [bottom_lim, top_lim]);
    else
        range_val = max_val - min_val_all;
        top_lim = max_val + 0.6 * range_val; % Headroom para leyenda
        bottom_lim = min_val_all - 0.1 * range_val; % Padding inferior
        ylim(ax, [bottom_lim, top_lim]);
    end

    % PADDING HORIZONTAL (X)
    % Ajustamos el margen lateral para que no toque los bordes izquierdo/derecho
    % Un margen de 0.8 unidades a cada lado suele ser ideal para grupos de barras
    xlim(ax, [0.2, num_active_groups + 0.8]);

    % --- 5. ESTILO DE EJES SIAM ---
    box(ax, 'on');
    grid(ax, 'on');
    ax.GridAlpha = 0.08;
    ax.Layer = 'top';

    set(ax, 'FontSize', 8, 'LineWidth', 0.6, 'TickLabelInterpreter', 'latex');

    ylabel(ax, y_label, 'Interpreter', 'latex', 'FontSize', 9);
    xlabel(ax, 'Problem', 'Interpreter', 'latex', 'FontSize', 9);
    % title(ax, sprintf('\\textbf{%s}', title_label), 'Interpreter', 'latex', 'FontSize', 10);

    xticks(ax, x_pos);
    xticklabels(ax, arrayfun(@(id) sprintf('%d', id), active_ids, 'UniformOutput', false));

    % --- 6. LEYENDA ---
    legend(ax, handles, legends, 'Interpreter', 'latex', ...
           'FontSize', 8, 'Location', 'northeast', 'Box', 'off');

    % --- 7. EXPORTACION TECNICA ---
    set(fig, 'Renderer', 'painters');
    clean_title = regexprep(title_label, '[^a-zA-Z0-9]', '_');
    export_path = fullfile(save_folder, sprintf('bar_plot_%s.eps', clean_title));
    print(fig, '-depsc2', '-r600', export_path);

    hold(ax, 'off');
end
