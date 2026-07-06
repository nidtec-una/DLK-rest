function plot_3D(data_cell, p)

    r_sphere = 20;
    % ============================
    % FIGURE – SIAM Proportions
    % ============================
    fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 17 8]);
    tlo = tiledlayout(1, length(data_cell), 'Padding', 'tight', 'TileSpacing', 'compact');

    % ============================
    % Determine global limits for consistency
    % ============================
    all_pmax = cellfun(@(x) max(x.p), data_cell);
    all_qmax = cellfun(@(x) max(x.q), data_cell);
    all_tmax = cellfun(@(x) max(x.t), data_cell);
    max_axes = max([all_pmax, all_qmax, all_tmax]) + 1;

    % Unique gradient
    g = gray(256);
    g = flipud(g(50:end - 30, :));
    colormap(fig, g);

    for i = 1:length(data_cell)
        ax = nexttile(tlo);
        hold(ax, 'on');
        T = data_cell{i};

        % ============================
        % Scatter 3D
        % ============================
        sc = scatter3(ax, T.p, T.q, T.t, r_sphere, T.Metric, ...
                      'o', 'filled', 'MarkerEdgeColor', [0.2 0.2 0.2], 'LineWidth', 0.4);
        sc.MarkerFaceAlpha = 0.8;
        sc.MarkerEdgeAlpha = 0.4;

        % ============================
        % Optimal point
        % ============================
        [m_val, idx] = max(T.Metric);
        p_opt = T.p(idx);
        q_opt = T.q(idx);
        t_opt = T.t(idx);
        plot3(ax, p_opt, q_opt, t_opt, 'ko', 'MarkerSize', 8, 'LineWidth', 1.2);

        % ============================
        % Bezier arrow
        % ============================
        p0 = [max_axes * 0.25, max_axes * 0.25, max_axes * 1.1]; % arrow origin
        dir = [p_opt, q_opt, t_opt] - p0;
        dir = dir / norm(dir);
        r_marker = 0.35 * mean([range(T.p), range(T.q), range(T.t)]) / 10;
        p3 = [p_opt, q_opt, t_opt] - r_marker * dir;
        p1 = p0 + 0.33 * (p3 - p0);
        p2 = p0 + 0.66 * (p3 - p0);
        s = linspace(0, 1, 80)';
        B = (1 - s).^3 .* p0 + 3 * (1 - s).^2 .* s .* p1 + 3 * (1 - s) .* s.^2 .* p2 + s.^3 .* p3;
        plot3(ax, B(:, 1), B(:, 2), B(:, 3), '--', 'LineWidth', 0.8, 'Color', [0.5 0.5 0.5]);

        % Arrowhead
        v = [p_opt, q_opt, t_opt] - B(end, :);
        v = v / norm(v);
        quiver3(ax, B(end, 1), B(end, 2), B(end, 3), 0.2 * v(1), 0.2 * v(2), 0.2 * v(3), ...
                0, 'LineWidth', 0.9, 'MaxHeadSize', 1.2, 'Color', [0.5 0.5 0.5]);

        % ============================
        % Text
        % ============================
        txt_str = sprintf('$\\theta_{max}(%d,%d,%d) = %.4f$', p_opt, q_opt, t_opt, m_val);
        text(ax, p0(1), p0(2), p0(3), txt_str, 'Interpreter', 'latex', ...
             'FontSize', 7.5, 'BackgroundColor', 'w', 'EdgeColor', [0.8 0.8 0.8], 'Margin', 2);

        % ============================
        % SIAM axes style
        % ============================
        view(ax, 35, 25);
        pbaspect(ax, [1 1 1]);
        grid(ax, 'on');
        ax.GridAlpha = 0.08;
        set(ax, 'FontSize', 8, 'LineWidth', 0.5, 'TickLabelInterpreter', 'latex', 'Box', 'on');
        xlabel(ax, '$p$', 'Interpreter', 'latex');
        ylabel(ax, '$q$', 'Interpreter', 'latex');
        zlabel(ax, '$t$', 'Interpreter', 'latex');
        title(ax, sprintf('$\\alpha = 10^{%d}$', round(log10(p.alphas(i)))), ...
              'Interpreter', 'latex', 'FontSize', 9);

        % ============================
        % Uniform limits
        % ============================
        xlim(ax, [0 max_axes]);
        ylim(ax, [0 max_axes]);
        zlim(ax, [0 max_axes]);
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
    % Global title (manual, on the figure)
    % ============================

    % axes(fig, 'Position',[0 0 1 1], 'Visible','off');  % Invisible axes for text
    % text(0.5, 0.9, sprintf('$\\theta(p,q,t)$ via DMDc-LQR-%s - %s', p.method, p.title), ...
    %    'HorizontalAlignment','center', 'Interpreter','latex', 'FontSize',10);

    % ============================
    % Export
    % ============================
    set(fig, 'Renderer', 'painters');
    print(fig, '-depsc2', '-r600', fullfile(p.path_results, sprintf('field_%s_plot_%s.eps', p.title, p.method)));
end
