function u = updatebyLQReco(Hbar, Bbar, xhat, m_last, alpha)
persistent T Hc Bc Kc_cached H_prev B_prev

try
    % --- Evitar recomputar si el sistema no cambió mucho ---
    if isempty(H_prev) || norm(H_prev - Hbar, 'fro') > 1e-6 || norm(B_prev - Bbar, 'fro') > 1e-6
        H_prev = Hbar; B_prev = Bbar;

        % --- Subespacio controlable reducido (QR en lugar de SVD) ---
        Co = ctrb(Hbar, Bbar);
        r = rank(Co);

        if r < size(Hbar,1)
            [Q,~,~] = svd(Co(:,1:r), 'econ');   % o usa qr(Co(:,1:r),0)
            T = Q; 
            Hc = T' * Hbar * T;
            Bc = T' * Bbar;

            Qc = eye(r);
            Rc = eye(size(Bc,2));
            [Kc_cached,~,~] = lqr(Hc, Bc, Qc, Rc);
        else
            Q = (1-alpha)*eye(size(Hbar));
            R = (alpha)*eye(size(Bbar, 2));
            [Kc_cached,~,~] = lqr(Hbar, Bbar, Q, R);
            T = eye(size(Hbar));  % No reducción
        end
    end

    % --- Aplicación del control ---
    xhatc = T' * xhat;
    u_opt = -Kc_cached * xhatc;
    u = floor(norm(u_opt));

catch
    warning('Fallo en LQR económico, se usa m anterior.');
    u = m_last;
end
