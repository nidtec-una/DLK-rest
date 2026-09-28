function run_tests(test_name)
    %
    %   run tests with code coverage
    %
    %   USAGE::
    %
    %       run_tests()                     % Runs all tests in tests/
    %       run_tests('test_performance')   % Runs specific test
    %       run_tests('test_collection.m')  % Runs specific test file
    %

    tic;

    matlab_dir = fileparts(mfilename('fullpath'));
    cd(matlab_dir);

    fprintf('\nHome is %s\n', getenv('HOME'));

    folder_to_cover = fullfile(matlab_dir, 'src');

    addpath(genpath(folder_to_cover));
    addpath(genpath(fullfile(matlab_dir, 'data')));
    addpath(genpath(fullfile(matlab_dir, 'tests')));

    % Determine target test file or folder
    if nargin < 1 || isempty(test_name)
        target_path = fullfile(matlab_dir, 'tests');
    else
        [~, name, ext] = fileparts(test_name);
        if isempty(ext)
            ext = '.m';
        end
        target_path = fullfile(matlab_dir, 'tests', [name, ext]);
        if ~exist(target_path, 'file') && exist(test_name, 'file')
            target_path = test_name;
        end
    end

    % Suppress warnings for nonexistent paths saved in MATLAB preferences
    warning('off', 'MATLAB:dispatcher:pathWarning');

    % Search for MOxUnit in local or neighboring directories
    moxunit_candidates = { ...
        fullfile(matlab_dir, 'MOxUnit'), ...
        fullfile(matlab_dir, '..', 'MOxUnit'), ...
        fullfile(matlab_dir, '..', '1.1. Repo Mejorado Primer Articulo', 'DLK-rest', 'MOxUnit'), ...
        fullfile(matlab_dir, '..', 'KrySBAS-MATLAB', 'MOxUnit') ...
    };

    for k = 1:length(moxunit_candidates)
        if exist(moxunit_candidates{k}, 'dir')
            addpath(genpath(moxunit_candidates{k}));
        end
    end

    % Initialize MOcov if available locally
    if exist(fullfile(matlab_dir, 'MOcov'), 'dir')
        addpath(genpath(fullfile(matlab_dir, 'MOcov')));
    end

    if exist('moxunit_runtests', 'file') == 0
        error('MOxUnit is not installed or not in the MATLAB path.');
    end

    if ispc
        success = moxunit_runtests(target_path, '-verbose');

    else
        success = moxunit_runtests(target_path, ...
                                   '-verbose', ...
                                   '-recursive', ...
                                   '-with_coverage', ...
                                   '-cover', ...
                                   folder_to_cover, ...
                                   '-cover_xml_file', ...
                                   'coverage.xml', ...
                                   '-cover_html_dir', ...
                                   fullfile(pwd, 'coverage_html') ...
                                  );
    end

    fileID = fopen('test_report.log', 'w');
    if success
        fprintf(fileID, '0');
    else
        fprintf(fileID, '1');
    end
    fclose(fileID);

    toc;

end
