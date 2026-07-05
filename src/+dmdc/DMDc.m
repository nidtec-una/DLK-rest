function [X_aprox_p,Xhat_p, Hbar, Bbar, Uhat, Htil, Btil] = DMDc(StateData, m, InputData, tau, q, t)
    [n, p_plus_p0] = size(StateData);
    p = p_plus_p0 - 1;

    % Minimizacion de ||R' - G*Omega|| con SVD truncado
    X = StateData(:,1:end-1);
    Xp = StateData(:, 2:end);
    Up = InputData(:, 1:end-1);
    Omega = [X; Up];

    [U, Sig, V] = svd(Omega, 'econ');
    fprintf('Rango truncado q = %d\n', q);
    Util = U(:, 1:q);
    Sigtil = Sig(1:q, 1:q);
    Vtil = V(:, 1:q);

    [U, ~, ~] = svd(StateData, 'econ');
    fprintf('Rango truncado para estados t = %d\n', t);
    Uhat = U(:, 1:t);

    U_1 = Util(1:n, :);
    U_2 = Util(n+1:n+1, :);
    Htil = Xp * Vtil * inv(Sigtil) * U_1';
    Btil = Xp * Vtil * inv(Sigtil) * U_2';
    Hbar = Uhat' *Htil* Uhat;
    Bbar = Uhat' * Btil;

    u = m;

    % valores de p
    Xhat_p(:, 1) = Uhat' * StateData(:, 1);
    for i = 1:p
        Xhat_p(:, i+1) = Hbar * Xhat_p(:, i) + Bbar * u;
    end

    try
        X_aprox_p = Uhat * Xhat_p;
    catch
        disp('Error en la multiplicación de matrices')
    end
    end