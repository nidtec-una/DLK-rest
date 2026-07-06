function T = calculate_data(A, b, p)
    % 1. Generate grid and filter valid points
    [P_grid, Q_grid, T_grid] = meshgrid(1:10, 1:10, 1:10);
    mask = (Q_grid <= P_grid) & (T_grid <= P_grid + 1);
    P = P_grid(mask);
    Q = Q_grid(mask);
    T_vec = T_grid(mask);

    num_pts = length(P);
    C_mean = zeros(num_pts, 1);
    C_rmssd = zeros(num_pts, 1);
    valid = false(num_pts, 1);

    % 2. Initialize progress bar with current method
    msg = sprintf('Calculating %s (alpha = %.0e)', p.method, p.current_alpha);
    h = waitbar(0, msg);

    for k = 1:num_pts
        try
            % DYNAMIC METHOD SELECTION
            if strcmpi(p.method, 'LGMRES')
                % [xf4, R4, totiter4, m_values4, p_values4, error_per_snp4] = snapshot_model.r.LGMRES_DMDc_LQR_fixp(A, b, x0, m0-l, l, tol, cycle_max,alpha);
                [~, ~, ~, m, ~, ~] = snapshot_model.r.LGMRES_DMDc_LQR_fixp( ...
                                                                      A, b, p.x0, p.m0 - p.l, p.l, p.tol, p.cycles, p.current_alpha, P(k), Q(k), T_vec(k));
            else % Default to GMRES
                [~, ~, ~, m, ~, ~] = snapshot_model.r.GMRES_DMDc_LQR_fixp( ...
                                                                     A, b, p.x0, p.m0, p.tol, p.cycles, p.current_alpha, P(k), Q(k), T_vec(k));
            end

            if length(m) > 2
                C_mean(k) = mean(m) / p.m_max;
                C_rmssd(k) = sqrt(mean(diff(m).^2)) / mean(m);
                valid(k) = true;
            end
        catch
            % Convergence error or interruption
        end

        % 3. Update progress bar
        if mod(k, 5) == 0 || k == num_pts
            waitbar(k / num_pts, h, sprintf('%s Alpha: %.0e | %d/%d', ...
                                          p.method, p.current_alpha, k, num_pts));
        end

        % Optional: Abort if waitbar is closed
        if ~ishandle(h)
            break
        end
    end

    % 4. Close bar and package results
    if ishandle(h)
        close(h);
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
