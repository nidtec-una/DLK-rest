function [xf, R, tot_cycles, m_values, p_values, erro_per_snp] = LGMRES_DMDc_LQR_fixp(A, b, x0, m0, l,tol, cycleMax,alpha_LQR,p,q,t)

%% Valores por defecto
if nargin < 11, t = 6; end
if nargin < 10, q = 5; end
if nargin < 9, p = 5; end
if nargin < 8, alpha_LQR = 1e-4; end
if nargin < 7, cycleMax = eps; end
if nargin < 6, tol = 1e-12; end
if nargin < 5, l = 3; end 
if nargin < 4, m0 = 27; end
if nargin < 3, x0 = zeros(size(b)); end
if nargin < 2, error('Se requieren al menos A y b'); end

%% Inicialización
tot_cycles = 0;
R = [];
m = m0;
m_values = [];
p_values = [];
xf = x0;

erro_per_snp = [];
j = 1; % contador de snapshots

while tot_cycles <= cycleMax && (isempty(R) || R(end)>= tol)
    
    m_values = [m_values, m];
    disp(['Snapshot ', num2str(j),', m = ', num2str(m)])
    
    %% Construir snapshot con LGMRES
    [MatR, xf, normRg, cycles_gmres] = sample_gmres.LGMRES(A, b, xf, m, l, tol, p);
    tot_cycles = tot_cycles + cycles_gmres;
    p_values = [p_values, cycles_gmres];

    if cycles_gmres < p
        R = [R, normRg];
        break;
    end

    normRg = normRg(1:end-1);
    R = [R, normRg];
    
    InputData = m * ones(1, size(MatR, 2));
    
    %% Identificación de sistema con DMDc
    tau = 0; % sin predicción
    [X_aprox_p, Xhat_p, Hbar, Bbar, ~, ~, ~] = dmdc.DMDc(MatR, m, InputData,tau, q,t);
    
    %% Calculo de error relativo por snapshot
    erro_snp = norm(MatR - X_aprox_p, 'fro') / norm(MatR, 'fro');
    erro_per_snp = [erro_per_snp, erro_snp];

    %% Calcular error por columna
    % Debe retirarse las utimas columnas para evitar duplicados
    %MatR_trim = MatR(:, 1:end-1);
    %X_aprox_trim = X_aprox_p(:, 1:end-1);
    %error_cols = vecnorm(MatR_trim - X_aprox_trim) ./ vecnorm(MatR_trim);
    %error = [error, error_cols];

    
    %% Actualización de m mediante LQR
    dim_sp = control_law.updatebyLQR(Hbar, Bbar, Xhat_p(:,end), m, alpha_LQR);
    % Acotado de m
    m_max = floor(size(A,2)*0.1);
    m = min(dim_sp-l, m_max-l);
    m_min = 5-l;
    m = max(m, m_min);

    j = j + 1; % incrementar contador de snapshots
    
end
end
