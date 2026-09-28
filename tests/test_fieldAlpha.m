function test_suite = test_fieldAlpha %#ok<*STOUT>
    try
        test_functions = localfunctions(); %#ok<*NASGU>
    catch
    end
    initTestSuite;
end

function test_fieldAlpha_matrix()
    % main_test_fieldAlpha.m

    % Configuration
    title = 'sherman5'; % Problem name
    curr_dir = fileparts(mfilename('fullpath'));
    mat_path = fullfile(curr_dir, '..', 'data', 'mat_collection_test_full', title);
    [A, b] = matrix_utils.load_matrices(mat_path);
    p.m0 = 30;
    p.l = 2;
    p.m_max = floor(size(A, 2) * 0.1);
    p.tol = 1e-9;
    p.x0 = zeros(size(b, 1), 1);
    p.cycles = 400;
    p.method = 'GMRES';
    p.alphas = [1e-4, 1e-5, 1e-6];
    p.title = title;
    p.folder = fullfile(curr_dir, '..', 'results', 'test_field');
    p.csv_folder = fullfile(p.folder, 'csv');

    if ~exist(p.folder, 'dir')
        mkdir(p.folder);
    end
    if ~exist(p.csv_folder, 'dir')
        mkdir(p.csv_folder);
    end
    results = cell(1, length(p.alphas));

    for i = 1:length(p.alphas)
        p.current_alpha = p.alphas(i);
        csv_file = fullfile(p.csv_folder, sprintf('res_%s_%s_a%.0e.csv', p.title, p.method, p.current_alpha));

        if exist(csv_file, 'file')
            fprintf(1, '\n[Alpha %d/%d] Loading cached: alpha %.0e (%s)...\n', i, length(p.alphas), p.current_alpha, p.title);
            drawnow;
            results{i} = readtable(csv_file);
        else
            fprintf(1, '\n[Alpha %d/%d] Calculating parameter sweep: alpha %.0e (%s)...\n', i, length(p.alphas), p.current_alpha, p.title);
            drawnow;
            res_table = utils.calculate_data(A, b, p);
            writetable(res_table, csv_file);
            results{i} = res_table;
        end
    end

    plots.plot_2D(results, p);
    assertTrue(length(results) == 3);
end
