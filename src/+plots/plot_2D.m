function plot_2D(data_cell, p)

    r_marker = 30;
    % ============================
    % FIGURE – SIAM Proportions (2D)
    % ============================
    fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 17 6.5]);
    tlo = tiledlayout(1, length(data_cell), 'Padding', 'tight', 'TileSpacing', 'compact');

    % Determine global limits for consistency
    all_pmax = cellfun(@(x) max(x.p), data_cell);
    all_tmax = cellfun(@(x) max(x.t), data_cell);
    max_axes = max([all_pmax, all_tmax]) + 1;

    % Unique gradient
    g = gray(256);
    g = flipud(g(50:end - 30, :));
    colormap(fig, g);

    for i = 1:length(data_cell)
        ax = nexttile(tlo);
        hold(ax, 'on');
        T = data_cell{i};

        % ============================
        % Scatter 2D (p vs t)
        % ============================
        sc = scatter(ax, T.p, T.t, r_marker, T.Metric, ...
                     'o', 'filled', 'MarkerEdgeColor', [0.2 0.2 0.2], 'LineWidth', 0.4);
        sc.MarkerFaceAlpha = 0.85;
        sc.MarkerEdgeAlpha = 0.5;

        % ============================
        % Optimal point
        % ============================
        [m_val, idx] = max(T.Metric);
        p_opt = T.p(idx);
        t_opt = T.t(idx);
        plot(ax, p_opt, t_opt, 'ko', 'MarkerSize', 8, 'LineWidth', 1.4, 'MarkerFaceColor', 'none');

        % Text annotation for max metric
        txt_str = sprintf('$\\theta_{max}(%d,%d) = %.4f$', p_opt, t_opt, m_val);
        text(ax, p_opt, t_opt + 0.8, txt_str, 'Interpreter', 'latex', ...
             'FontSize', 7.5, 'BackgroundColor', 'w', 'EdgeColor', [0.8 0.8 0.8], ...
             'Margin', 2, 'HorizontalAlignment', 'center');

        % ============================
        % SIAM axes style
        % ============================
        grid(ax, 'on');
        ax.GridAlpha = 0.1;
        set(ax, 'FontSize', 8, 'LineWidth', 0.5, 'TickLabelInterpreter', 'latex', 'Box', 'on');
        xlabel(ax, '$p = q$', 'Interpreter', 'latex');
        ylabel(ax, '$t$', 'Interpreter', 'latex');
        title(ax, sprintf('$\\alpha = 10^{%d}$', round(log10(p.alphas(i)))), ...
              'Interpreter', 'latex', 'FontSize', 9);

        xlim(ax, [0 max_axes]);
        ylim(ax, [0 max_axes]);
    end

    % ============================
    % Single colorbar
    % ============================
    cb = colorbar;
    cb.Layout.Tile = 'east';
    cb.Label.String = 'Metric Value ($\theta$)';
    cb.Label.Interpreter = 'latex';
    cb.TickLabelInterpreter = 'latex';
    cb.FontSize = 8;

    % ============================
    % Export
    % ============================
    eps_folder = fullfile(p.folder, 'eps');
    png_folder = fullfile(p.folder, 'png');
    if ~exist(eps_folder, 'dir')
        mkdir(eps_folder);
    end
    if ~exist(png_folder, 'dir')
        mkdir(png_folder);
    end
    if ishandle(fig)
        saveas(fig, fullfile(eps_folder, sprintf('field_%s_plot_2D.eps', p.title)), 'epsc');
        saveas(fig, fullfile(png_folder, sprintf('field_%s_plot_2D.png', p.title)));
        close(fig);
    end
end
