% --- Ley de control con Riccati con lqr + svd---
function u = updatebyLQR(Hbar, Bbar, xhat, m_last, alpha)

    Co = ctrb(Hbar, Bbar);
    r = rank(Co);

    if r < size(Hbar, 1)
        [U, ~, ~] = svd(Co, 'econ');
        T = orth(U(:, 1:r));      % Base ortonormal del subespacio controlable

        Hc = T' * Hbar * T;      % Sistema reducido en coordenadas controlables
        Bc = T' * Bbar;

        try
            Qc = eye(r);
            Rc = eye(size(Bc, 2));
            [Kc, ~, ~] = lqr(Hc, Bc, Qc, Rc);
            xhatc = T' * xhat; % Estado en coordenadas controlables
            xhatc = xhatc / norm(xhatc);
            u_opt = -Kc * xhatc;
            u = floor(norm(u_opt));
        catch
            warning('LQR también falló en el subespacio controlable. Se usará la m anterior.');
            u = m_last;
            return
        end
    else
        try
            Q = (1 - 0) * eye(size(Hbar));  %% modificado por jccf
            R = (alpha) * eye(size(Bbar, 2)); %% modificado por jccf
            [K, ~, ~] = lqr(Hbar, Bbar, Q, R);
            xhat = xhat / norm(xhat);
            u_opt = -K * xhat;
            u = floor(norm(u_opt));
        catch
            warning('El sistema es completamente controlable, pero LQR falló. Se usará la m anterior.');
            u = m_last;
            return
        end
    end
