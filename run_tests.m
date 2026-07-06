function run_tests()
    %
    %   run tests with code coverage
    %
    %   USAGE::
    %
    %       run_tests()
    %

    tic;

    matlab_dir = fileparts(mfilename('fullpath'));
    cd(matlab_dir);

    fprintf('\nHome is %s\n', getenv('HOME'));

    folder_to_cover = fullfile(matlab_dir, 'src');

    test_folder = fullfile(matlab_dir, 'tests');

    addpath(genpath(folder_to_cover));
    addpath(genpath(fullfile(matlab_dir, 'data')));

    % Initialize MOxUnit if available locally (useful for GitHub Actions)
    if exist(fullfile(matlab_dir, 'MOxUnit'), 'dir')
        addpath(genpath(fullfile(matlab_dir, 'MOxUnit')));
    end

    % Initialize MOcov if available locally
    if exist(fullfile(matlab_dir, 'MOcov'), 'dir')
        addpath(genpath(fullfile(matlab_dir, 'MOcov')));
    end

    if exist('moxunit_runtests', 'file') == 0
        error('MOxUnit is not installed or not in the MATLAB path.');
    end

    if ispc
        success = moxunit_runtests(test_folder, '-verbose');

    else
        success = moxunit_runtests(test_folder, ...
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
