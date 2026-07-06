% --- Ley de control con PD ---
function u = updatebyPD(Rdmd, m)
    u = m - floor(-3 * Rdmd(end) / Rdmd(end - 1)) + floor(9 * (Rdmd(end) - Rdmd(end - 2)) / (2 * Rdmd(end - 1)));
end
