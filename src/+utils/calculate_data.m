function T = calculate_data(A, b, p)
    % 1. Generate 2D grid (p = q) and filter valid points
    [P_grid, T_grid] = meshgrid(1:10, 1:10);
    mask = (T_grid <= P_grid + 1);
    P = P_grid(mask);
    Q = P; % p = q
    T_vec = T_grid(mask);

    num_pts = length(P);
    C_mean = zeros(num_pts, 1);
    C_rmssd = zeros(num_pts, 1);
    valid = false(num_pts, 1);

    % 2. Initialize progress / loading indicator
    has_gui = usejava('awt');
    hWait = [];
    msg = sprintf('Calculating DMDc-LQR-GMRES (alpha = %.0e)', p.current_alpha);
    if has_gui
        try
            hWait = waitbar(0, msg, 'Name', 'Parameter Sweep Progress');
        catch
            hWait = [];
        end
    end
    fprintf('\n%s\n', msg);

    for k = 1:num_pts
        try
            [~, ~, ~, m, ~, ~] = snapshot_model.r.GMRES_DMDc_LQR_fixp( ...
                                                                 A, b, p.x0, p.m0, p.tol, p.cycles, p.current_alpha, P(k), Q(k), T_vec(k));

            if length(m) > 2
                C_mean(k) = mean(m) / p.m_max;
                C_rmssd(k) = sqrt(mean(diff(m).^2)) / mean(m);
                valid(k) = true;
            end
        catch
            % Convergence error or interruption
        end

        % 3. Update progress bar (Console + GUI waitbar)
        progress = k / num_pts;
        pct = round(progress * 100);
        num_bars = round(progress * 20);
        bar_str = [repmat('=', 1, num_bars), repmat('-', 1, 20 - num_bars)];

        if mod(k, 2) == 0 || k == num_pts
            fprintf(1, '\rProgress: [%s] %3d%% (%d/%d)', bar_str, pct, k, num_pts);
            drawnow;
            if ~isempty(hWait) && ishandle(hWait)
                waitbar(progress, hWait, sprintf('Alpha: %.0e | [%s] %d%%', p.current_alpha, bar_str, pct));
            end
        end
    end
    fprintf('\nDone.\n');

    % 4. Close GUI waitbar if open
    if ~isempty(hWait) && ishandle(hWait)
        close(hWait);
    end

    % Filter only valid data
    T = table(P(valid), Q(valid), T_vec(valid), C_mean(valid), C_rmssd(valid), ...
              'VariableNames', {'p', 'q', 't', 'mean', 'rmssd'});

    % 5. Metric calculations (Protection against empty table)
    if ~isempty(T)
        T.S_norm = 1 - (T.rmssd / max(T.rmssd));
        T.Centrality = max(0, 1 - abs((T.mean * p.m_max) - p.m_max / 2) / (p.m_max / 2));
        T.Metric = (T.S_norm + T.Centrality) / 2;
    end
end
