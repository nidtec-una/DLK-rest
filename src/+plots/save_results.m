% save_results guarda los resultados en un archivo
%   save_results(m, tol, pmin, pmax, x0, Rg, Rdmd, R_tau, title_label)
%   guarda los resultados en un archivo
%   m: orden del sistema
%   tol: tolerancia
%   pmin: minimo p
%   pmax: maximo p
%   x0: vector de inicialización
%   Rg: vector de normas de los errores DMDc
%   Rdmd: vector de normas de los errores DMDc con autocorrelación
%   R_tau: vector de normas de los errores DMDc con autocorrelación
%   title_label: etiqueta del archivo
function save_results(m, tol, pmin, pmax, x0, Rg, Rdmd, R_tau, title_label)
    if ~exist('test', 'dir')
        mkdir('test');
    end
    clean_title = regexprep(title_label, '[^a-zA-Z0-9]', '');
    fileID = fopen(fullfile('test', [clean_title, '.txt']), 'w');
    fprintf(fileID, 'm: %d\ntol: %e\npmin: %d\npmax: %d\n', m, tol, pmin, pmax);
    fprintf(fileID, 'x0: [%s]\n', num2str(x0'));
    fclose(fileID);
end
