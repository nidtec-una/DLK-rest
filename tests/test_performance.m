function test_suite = test_performance %#ok<*STOUT>
    try
        test_functions = localfunctions(); %#ok<*NASGU>
    catch
    end
    initTestSuite;
end

function test_performance_sherman5()
    % main_test_performance.m logic
    % --- Parameters ---
    curr_dir = fileparts(mfilename('fullpath'));
    mat_path = fullfile(curr_dir, '..', 'data', 'mat_collection_test_full', 'sherman5');
    [A, b] = matrix_utils.load_matrices(mat_path);
    title_label = 'sherman5';
    m0 = 30;
    m = 110;
    % [A,b] = matrix_utils.genMatrixZNn(50,100);
    % title_label = 'Zhong & Morgan';
    % m0 = 30;

    tol = 1e-9;
    x0 = zeros(size(b, 1), 1);
    normR0 = norm(b - A * x0);
    cycle_max = 400;
    alpha = 1e-4; %%% modificado por jccf

    fprintf(1, '\nRunning test_performance on matrix %s...\n', title_label);
    % --- GMRES-DMDc-LQR con p fijo ---
    fprintf(1, '[1/3] Calculating DMDc-LQR-GMRES...\n');
    drawnow;
    tic;
    [xf3, R3, totiter3, m_values3, p_values3, error_per_snp3] = snapshot_model.r.GMRES_DMDc_LQR_fixp(A, b, x0, m0, tol, cycle_max, alpha);
    t_dmdc_gmres = toc;

    % --- GMRES ---
    fprintf(1, '[2/3] Calculating standard GMRES...\n');
    drawnow;
    tic;
    [~, ~, Rgmres, cycles_gmres] = sample_gmres.GMRES(A, b, x0, m, tol, cycle_max);
    t_gmres = toc;

    % --- LGMRES ---
    fprintf(1, '[3/3] Calculating standard LGMRES...\n');
    drawnow;
    l = 3;
    m_lgmres = m - l;
    tic;
    [~, ~, Rlgmres, cycles_lgmres] = sample_gmres.LGMRES(A, b, x0, m_lgmres, l, tol, cycle_max);
    t_lgmres = toc;

    color_gmres   = [0.1, 0.4, 0.8];   % Azul
    color_lgmres  = [1.0, 0.5, 0.0];   % Naranja
    color_R3      = [0.3, 0.7, 0.5];   % Verde azulado

    R_list = {
              {Rgmres, ['GMRES(', num2str(m), '), cycles = ', num2str(cycles_gmres)], ':o', color_gmres}
              {Rlgmres, ['LGMRES(', num2str(m_lgmres), ',', num2str(l), '), cycles = ', num2str(cycles_lgmres)], '--s', color_lgmres}
              {R3, ['DMDc-LQR-GMRES, cycles = ', num2str(totiter3), ')'], '--^', color_R3}
             };

    % -- Methods comparison plot
    save_folder = fullfile(curr_dir, '..', 'results', 'test_performance');
    if ~exist(save_folder, 'dir')
        mkdir(save_folder);
    end
    plots.show_results_custom(R_list, ['Comparison of methods - ', title_label], tol, save_folder);
    m_limits = [2, floor(size(A, 2) * 0.1)];

    % -- DMDc-LQR-GMRES performance plot
    R_gmres_dmdc = {
                    {R3, error_per_snp3,     'DMDcR-LQR-GMRES', '--^', [0.4, 0.4, 0.4]}
                   };
    plots.show_results_with_m(R_gmres_dmdc, ['Performance of DMDc-LQR-GMRES($m$) - ', title_label], m_limits, m_values3, p_values3, save_folder);

    % -- Export CSV results --
    csv_folder = fullfile(save_folder, 'csv');
    if ~exist(csv_folder, 'dir')
        mkdir(csv_folder);
    end
    perf_data = {
        'Method', 'Cycles', 'Time_sec', 'Final_Residual';
        sprintf('GMRES(%d)', m), cycles_gmres, t_gmres, Rgmres(end);
        sprintf('LGMRES(%d,%d)', m_lgmres, l), cycles_lgmres, t_lgmres, Rlgmres(end);
        'DMDc-LQR-GMRES', totiter3, t_dmdc_gmres, R3(end)
    };
    writecell(perf_data, fullfile(csv_folder, [title_label, '_performance.csv']));

    % Asserts for MOxUnit
    assertTrue(totiter3 > 0);
end
