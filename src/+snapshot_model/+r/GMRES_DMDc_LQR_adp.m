function [xf, R, tot_cycles, m_values, p_values, error_per_snp, Disp] = GMRES_DMDc_LQR_adp(A, b, x0, m0, tol, cycleMax, p_min, p_max, threshSVD,alpha_LQR)
% GMRES_DMDc_LQR_adp3 - GMRES adaptativo con DMDc y actualización LQR
% Esta versión adapta p usando la diferencia inmediata en escala logarítmica:
% Disp_raw = | log10(normRg(end)) - log10(normRg(1)) |
% p se mapea invertido: baja dispersión -> p grande, alta dispersión -> p pequeño
%
% Firma:
% [xf, R, tot_cycles, m_values, p_values, Disp] = GMRES_DMDc_LQR_adp3(A, b, x0, m0, tol, cycleMax, p_min, p_max, threshSVD)

%% Valores por defecto
if nargin < 10, alpha_LQR = 1e-4; end
if nargin < 9, threshSVD = eps; end
if nargin < 8, p_max = 6; end
if nargin < 7, p_min = 3; end
if nargin < 6, cycleMax = 300; end
if nargin < 5, tol = eps; end
if nargin < 4, m0 = 30; end
if nargin < 3, x0 = zeros(size(b)); end
if nargin < 2, error('Se requieren al menos A y b'); end


%% Inicialización
tot_cycles = 0;
R = [];
normR0 = norm(b - A * x0);
m = m0;
m_values = [];
p_values = [];
xf = x0;

Disp = [];
error_per_snp = [];

alpha = 0.1;        % factor de suavizado EMA
Disp_smooth_prev = []; 

j = 1; % contador de snapshots/iteraciones

while tot_cycles <= cycleMax && (isempty(R) || R(end)/normR0 >= tol)

    m_values = [m_values, m];

    %% Selección de p usando la dispersión suavizada previa
    if j == 1
        p = p_min; % inicialización conservadora
    else
        factor = 1 - Disp_smooth_prev; % invertido: baja Disp -> p grande
        p = p_min + ceil((p_max - p_min) * factor);
        p = max(p_min, min(p_max, p)); % límites
    end

    disp(['Snapshot ', num2str(j), ': p = ', num2str(p), ', m = ', num2str(m)])

    %% Construir snapshot con GMRES
    [MatR, xf, normRg, cycles_gmres] = sample_gmres.GMRES(A, b, xf, m, tol, p);
    tot_cycles = tot_cycles + cycles_gmres;
    p_values = [p_values, cycles_gmres];
    
    % Condicion de parada si ocurre convergencia prematura
    if cycles_gmres < p
        R = [R, normRg];
        break;
    end

    normRg = normRg(1:end-1);
    R = [R, normRg];
    
    InputData = m * ones(1, size(MatR, 2));
    
    %% --- Dispersión logarítmica ---
    r1 = max(normRg(1), eps);
    rk = max(normRg(end), eps);
    Disp_raw = abs(log10(rk) - log10(r1)); % diferencia logarítmica
    % Normalizar para que quede en [0,1] aproximadamente
    Disp_norm = min(1, Disp_raw / 5);  % escala: 1 unidad log ≈ factor 10
    Disp = [Disp, Disp_norm];
    
    %% Suavizado exponencial (EMA)
    if j == 1
        Disp_smooth = Disp_norm;
    else
        Disp_smooth = alpha * Disp_norm + (1 - alpha) * Disp_smooth_prev;
    end
    Disp_smooth_prev = Disp_smooth;

    %% Identificación de sistema con DMDc
    tau = 0; % sin predicción
    [X_aprox_p, Xhat_p, Hbar, Bbar, ~, ~, ~] = dmdc.DMDc(MatR, m, InputData, tau, threshSVD);

    %% Calculo de error relativo por snapshot
    erro_snp = norm(MatR - X_aprox_p, 'fro') / norm(MatR, 'fro');
    error_per_snp = [error_per_snp, erro_snp];

    %% Actualización de m mediante LQR
    m = control_law.updatebyLQR(Hbar, Bbar, Xhat_p(:,end), m, alpha_LQR);
    % Acotado de m
    m_max = floor(size(A,2)*0.1);
    m = min(m, m_max);
    m_min = 2;
    m = max(m, m_min);

    j = j + 1; % incrementar contador de snapshots
end
end
