% main_test_fieldAlpha.m
clc; clear; close all;
addpath(genpath('funciones'));

% Configuración
title = 'add20'; % Nombre del problema
[A, b] = construir_matriz.cargar_matices(strcat('mat_collection_test_full/', title));
p.title = title;
p.m0 = 30;
p.l = 2;
p.m_max = floor(size(A,2)*0.1);
p.tol = 1e-9;
p.x0 = zeros(size(b,1),1);
p.cycles = 400;
p.method = 'LGMRES';
p.alphas = [1e-4, 1e-5, 1e-6];
p.folder = 'test_field';
p.path_results = fullfile(p.folder, 'eps');

if ~exist(p.folder, 'dir'), mkdir(p.folder); end
resultados = cell(1, length(p.alphas));

for i = 1:length(p.alphas)
    p.current_alpha = p.alphas(i);
    csv_file = fullfile(p.folder, sprintf('res_%s_%s_a%.0e.csv', p.title, p.method, p.current_alpha));
    
    if exist(csv_file, 'file')
        fprintf('Cargando: alpha %.0e\n', p.current_alpha);
        resultados{i} = readtable(csv_file);
    else
        fprintf('Calculando: alpha %.0e\n', p.current_alpha);
        res_table = utils.calcular_datos(A, b, p);
        writetable(res_table, csv_file);
        resultados{i} = res_table;
    end
end

graficos.graficar_3D(resultados, p);