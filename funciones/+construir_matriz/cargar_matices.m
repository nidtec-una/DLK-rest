% cargar_matrices carga una matriz de la carpeta mat_collection
%   [A, b] = cargar_matrices(ruta_archivo)
%   carga una matriz de la carpeta mat_collection
%   ruta_archivo: ruta de la matriz
%   A: matriz   
%   b: vector de las entradas
function [A, b] = cargar_matrices(ruta_archivo)
    load(ruta_archivo);
    % if load b fail load a vector of ones
    try
        b = Problem.b;
    catch
        b = ones(size(Problem.A,1),1);
    end
    A = Problem.A;
    end
    