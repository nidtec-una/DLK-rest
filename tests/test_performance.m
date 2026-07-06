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

    % --- GMRES-DMDc-LQR con p fijo ---
    tic;
    % p = 5;
    % q = 5;
    % t = 6;
    [xf3, R3, totiter3, m_values3, p_values3, error_per_snp3] = snapshot_model.r.GMRES_DMDc_LQR_fixp(A, b, x0, m0, tol, cycle_max, alpha);
    t_dmdc_gmres = toc;

    % --- LGMRES-DMDcR-LQR con p fijo ---
    tic;
    l = 2;
    % p = 10;
    % q = 2;
    % t = 7;
    [xf4, R4, totiter4, m_values4, p_values4, error_per_snp4] = snapshot_model.r.LGMRES_DMDc_LQR_fixp(A, b, x0, m0 - l, l, tol, cycle_max, alpha);
    t_dmdc_lgmres = toc;

    % --- GMRES ---
    tic;
    [~, ~, Rgmres, cycles_gmres] = sample_gmres.GMRES(A, b, x0, m, tol, cycle_max);
    t_gmres = toc;

    % --- LGMRES ---
    l = 3;
    m_lgmres = m - l;
    tic;
    [~, ~, Rlgmres, cycles_lgmres] = sample_gmres.LGMRES(A, b, x0, m_lgmres, l, tol, cycle_max);
    t_lgmres = toc;

    color_gmres   = [0.1, 0.4, 0.8];   % Azul
    color_lgmres  = [1.0, 0.5, 0.0];   % Naranja
    color_R3      = [0.3, 0.7, 0.5];   % Verde azulado
    color_R4      = [0.6, 0.2, 0.7];  % Morado

    R_list = {
              {Rgmres, ['GMRES(', num2str(m), '), cycles = ', num2str(cycles_gmres)], ':o', color_gmres}
              {Rlgmres, ['LGMRES(', num2str(m_lgmres), ',', num2str(l), '), cycles = ', num2str(cycles_lgmres)], '--s', color_lgmres}
              {R3, ['DMDc-LQR-GMRES, cycles = ', num2str(totiter3), ')'], '--^', color_R3}
              {R4, ['DMDc-LQR-LGMRES, cycles = ', num2str(totiter4), ')'], '-d', color_R4}
             };

    % -- Methods comparison plot
    save_folder = fullfile(curr_dir, 'results', 'test_performance');
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

    % -- DMDc-LQR-LGMRES performance plot
    R_gmres_dmdc = {
                    {R4, error_per_snp4,     'DMDcR-LQR-LGMRES', '-d', [0.6, 0.6, 0.6]}
                   };
    % plots.show_results_with_m(R_gmres_dmdc, ['Performance of DMDc-LQR-LGMRES($m,\ell$) - ', title_label], m_limits, m_values4, p_values4, save_folder);

    % Asserts for MOxUnit
    assertTrue(totiter3 > 0);
    assertTrue(totiter4 > 0);
end
