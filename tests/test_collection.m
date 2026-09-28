function test_suite = test_collection %#ok<*STOUT>
    try
        test_functions = localfunctions(); %#ok<*NASGU>
    catch
    end
    initTestSuite;
end

function test_collection_matrices()
    % main_test_collection.m

    % --- Folders and Files ---
    curr_dir = fileparts(mfilename('fullpath'));
    mat_folder = fullfile(curr_dir, '..', 'data', 'mat_collection_test_full');
    save_folder = fullfile(curr_dir, '..', 'results', 'test_collection');
    csv_folder = fullfile(save_folder, 'csv');
    results_file = fullfile(csv_folder, 'results_avg_full.csv');

    if ~exist(save_folder, 'dir')
        mkdir(save_folder);
    end

    % Define general parameters
    mat_files = dir(fullfile(mat_folder, '*.mat'));
    num_matrices = length(mat_files);

    %% --- CONTROL LOGIC: RUN OR LOAD? ---
    if exist(results_file, 'file') && size(readcell(results_file), 1) > 1
        fprintf('CSV files already exist. Loading data to plot directly...\n');
        results = readcell(results_file);

    else
        fprintf('No previous results found. Starting calculations...\n');

        % --- Simulation parameters ---
        tol = 1e-9;
        max_cycles = 400;
        m0 = 30;
        l = 2;
        m_list = repmat(30, 1, num_matrices);
        num_runs = 1;

        % --- Initialize containers ---
        matrix_info = { 'Problem', 'Matrix', 'n', 'nnz', 'cond_A' };
        cycles_csv = { 'Problem', 'Matrix', 'Cycles_GMRES', 'Cycles_LGMRES', 'Cycles_GMRES_DMDc_LQR' };
        rel_time_csv = { 'Problem', 'Matrix', 'avg_time_GMRES', 'avg_time_LGMRES', 'avg_time_GMRES_DMDc_LQR' };
        results = {'Problem', 'Matrix', 'n', 'nnz', 'm_GMRES', 'm_l_GMRES', ...
                   'Cycles_GMRES', 'Cycles_LGMRES', 'Cycles_GMRES_DMDc_LQR', ...
                   'Avg_Time_GMRES', 'Avg_Time_LGMRES', 'Avg_Time_GMRES_DMDc_LQR'};

        problem_id = 1;
        % Each matrix has 3 method executions per run
        total_steps = num_matrices * num_runs * 3;
        current_step = 0;

        has_gui = usejava('awt');
        hWait = [];
        if has_gui
            try
                hWait = waitbar(0, 'Initializing experiment...', 'Name', 'Running Collection Experiments');
            catch
                hWait = [];
            end
        end

        for k = 1:num_matrices
            mat_name = mat_files(k).name;
            mat_name_str = strrep(erase(mat_name, '.mat'), '_', '\_');
            mat_path = fullfile(mat_folder, mat_name);

            [A, b] = matrix_utils.load_matrices(mat_path);
            n = size(A, 1);
            nnzA = nnz(A);
            x0 = zeros(n, 1);
            condA = condest(A);
            condA_str = sprintf('%.4e', condA);
            m = m_list(k);

            times_gmres = zeros(num_runs, 1);
            times_lgmres = zeros(num_runs, 1);
            times_gmres_dmdc = zeros(num_runs, 1);

            for run = 1:num_runs
                % Helper function to update loading progress
                update_wb = @(method) update_collection_progress(current_step, total_steps, k, num_matrices, mat_name, method, hWait);

                % 1. GMRES
                current_step = current_step + 1;
                update_wb('GMRES');
                tic;
                [~, ~, Rgmres, cycles_gmres] = sample_gmres.GMRES(A, b, x0, m, tol, max_cycles);
                times_gmres(run) = toc;
                last_R_gmres = Rgmres;

                % 2. LGMRES
                current_step = current_step + 1;
                update_wb('LGMRES');
                tic;
                [~, ~, Rlgmres, cycles_lgmres] = sample_gmres.LGMRES(A, b, x0, m - l, l, tol, max_cycles);
                times_lgmres(run) = toc;
                last_R_lgmres = Rlgmres;

                % 3. GMRES-DMDc
                current_step = current_step + 1;
                update_wb('GMRES-DMDc');
                tic;
                [~, R3, iter3, ~] = snapshot_model.r.GMRES_DMDc_LQR_fixp(A, b, x0, m0, tol, max_cycles);
                times_gmres_dmdc(run) = toc;
                last_R_gmres_dmdc = R3;
            end

            % Averages and convergence check
            avg_gmres = mean(times_gmres);
            avg_lgmres = mean(times_lgmres);
            avg_gmres_dmdc = mean(times_gmres_dmdc);

            cg = cycles_gmres;
            cl = cycles_lgmres;
            cgD = iter3;

            if isempty(last_R_gmres) || last_R_gmres(end) >= tol
                cg = NaN;
                avg_gmres = NaN;
            end
            if isempty(last_R_lgmres) || last_R_lgmres(end) >= tol
                cl = NaN;
                avg_lgmres = NaN;
            end
            if isempty(last_R_gmres_dmdc) || last_R_gmres_dmdc(end) >= tol
                cgD = NaN;
                avg_gmres_dmdc = NaN;
            end

            % Save in results structure
            results(problem_id + 1, :) = {problem_id, mat_name_str, n, nnzA, m, sprintf('%d,%d', m - l, l), ...
                                       cg, cl, cgD, avg_gmres, avg_lgmres, avg_gmres_dmdc};

            matrix_info(end + 1, :) = {problem_id, mat_name_str, n, nnzA, condA_str};

            % Robust cleanup function (avoids char conversion error)
            f_csv = @(val) string(val).replace("NaN", "NC");

            cycles_csv(end + 1, :) = {problem_id, mat_name_str, f_csv(cg), f_csv(cl), f_csv(cgD)};
            rel_time_csv(end + 1, :) = {problem_id, mat_name_str, f_csv(avg_gmres), f_csv(avg_lgmres), f_csv(avg_gmres_dmdc)};

            problem_id = problem_id + 1;
        end
        if ishandle(hWait)
            close(hWait);
        end

        % Export
        if ~exist(csv_folder, 'dir')
            mkdir(csv_folder);
        end
        writecell(matrix_info, fullfile(csv_folder, 'matrix_info.csv'));
        writecell(cycles_csv, fullfile(csv_folder, 'cycles_experiment.csv'));
        writecell(rel_time_csv, fullfile(csv_folder, 'avg_relativa_time_experiment.csv'));
        writecell(results, results_file);
    end

    %% --- DATA PROCESSING FOR PLOTS ---
    x_labels = results(2:end, 2);
    raw_data = cell2mat(cellfun(@(x) double(string(x).replace("NC", "NaN")), results(2:end, 7:12), 'UniformOutput', false));

    cycles_gmres_f       = raw_data(:, 1);
    cycles_lgmres_f      = raw_data(:, 2);
    cycles_gmres_dmdc_f  = raw_data(:, 3);

    avg_time_gmres_f       = raw_data(:, 4);
    avg_time_lgmres_f      = raw_data(:, 5);
    avg_time_gmres_dmdc_f  = raw_data(:, 6);

    %% --- PLOT CONFIGURATION ---
    gray_colors = [0.2 0.2 0.2; 0.5 0.5 0.5];

    time_list = {
                 {log10(avg_time_gmres_f ./ avg_time_lgmres_f),      'log(GMRES$(30)$/LGMRES$(27,3)$)',      gray_colors(1, :)}, ...
                 {log10(avg_time_gmres_f ./ avg_time_gmres_dmdc_f),  'log(GMRES$(30)$/DMDc-LQR-GMRES)',  gray_colors(2, :)}
                };

    cycle_list = {
                  {log10(cycles_gmres_f ./ cycles_lgmres_f),      'log(GMRES$(30)$/LGMRES$(27,3)$)',      gray_colors(1, :)}, ...
                  {log10(cycles_gmres_f ./ cycles_gmres_dmdc_f),  'log(GMRES$(30)$/DMDc-LQR-GMRES)',  gray_colors(2, :)}
                 };

    plots.show_bar_dual_custom(time_list, 'Average computation relative time', 'log(GMRES($30$)/Krylov Solver)', false, save_folder, x_labels);
    plots.show_bar_dual_custom(cycle_list, 'Number of cycles','log(GMRES($30$)/Krylov Solver)', false, save_folder, x_labels);

    fprintf('Process finished. Plots generated in: %s\n', save_folder);
    assertTrue(exist(save_folder, 'dir') == 7);
end

function update_collection_progress(current_step, total_steps, k, num_matrices, mat_name, method, hWait)
    progress = current_step / total_steps;
    pct = round(progress * 100);
    num_bars = round(progress * 20);
    bar_str = [repmat('=', 1, num_bars), repmat('-', 1, 20 - num_bars)];
    
    fprintf(1, '\rProgress: [%s] %3d%% | Matrix %d/%d (%s) [%s]', bar_str, pct, k, num_matrices, mat_name, method);
    drawnow;
    
    if ~isempty(hWait) && ishandle(hWait)
        waitbar(progress, hWait, sprintf('Matrix %d/%d: %s [%s] (%d%%)', k, num_matrices, mat_name, method, pct));
    end
end
