function test_suite = test_fieldAlpha %#ok<*STOUT>
    try
        test_functions = localfunctions(); %#ok<*NASGU>
    catch
    end
    initTestSuite;
end

function test_fieldAlpha_add20()
    % main_test_fieldAlpha.m

    % Configuration
    title = 'add20'; % Problem name
    curr_dir = fileparts(mfilename('fullpath'));
    mat_path = fullfile(curr_dir, '..', 'data', 'mat_collection_test_full', title);
    [A, b] = matrix_utils.load_matrices(mat_path);
    p.m0 = 30;
    p.l = 2;
    p.m_max = floor(size(A, 2) * 0.1);
    p.tol = 1e-9;
    p.x0 = zeros(size(b, 1), 1);
    p.cycles = 400;
    p.method = 'LGMRES';
    p.alphas = [1e-4, 1e-5, 1e-6];
    p.title = title;
    p.folder = fullfile(curr_dir, 'results', 'test_field');
    p.path_results = fullfile(p.folder, 'eps');

    if ~exist(p.folder, 'dir')
        mkdir(p.folder);
    end
    results = cell(1, length(p.alphas));

    for i = 1:length(p.alphas)
        p.current_alpha = p.alphas(i);
        csv_file = fullfile(p.folder, sprintf('res_%s_%s_a%.0e.csv', p.title, p.method, p.current_alpha));

        if exist(csv_file, 'file')
            fprintf('Loading: alpha %.0e\n', p.current_alpha);
            results{i} = readtable(csv_file);
        else
            fprintf('Calculating: alpha %.0e\n', p.current_alpha);
            res_table = utils.calculate_data(A, b, p);
            writetable(res_table, csv_file);
            results{i} = res_table;
        end
    end

    % plots.plot_3D(results, p);
    assertTrue(length(results) == 3);
end
