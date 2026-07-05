% load_matrices loads a matrix and vector from the specified path
%
%   [A, b] = load_matrices(file_path)
%
%   Inputs:
%       file_path: path to the .mat file
%   Outputs:
%       A: sparse matrix   
%       b: right-hand side vector
function [A, b] = load_matrices(file_path)
    load(file_path);
    % If loading b fails, generate a vector of ones
    try
        b = Problem.b;
    catch
        b = ones(size(Problem.A,1),1);
    end
    A = Problem.A;
    end
    