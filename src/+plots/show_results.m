function show_results(Rg, Rginic, Rdmd_f, Rdmd_a, R_tau, m, p, title_label)
    figure('Units', 'centimeters', 'Position', [2, 2, 14, 10]);

    semilogy(Rg, '-^', 'LineWidth', 2);
    hold on;
    semilogy(Rdmd_f, '--o', 'LineWidth', 1);
    semilogy(Rdmd_a, '--s', 'LineWidth', 1);

    title(title_label, 'Interpreter', 'latex', 'FontSize', 12);

    xlabel('Number of cycles', 'Interpreter', 'latex', 'FontSize', 10);
    ylabel('Residual norm', 'Interpreter', 'latex', 'FontSize', 10);

    legend({['GMRES(', num2str(m), ')'], ['PD-DMDc(p=', num2str(p), ')'], ['Riccati-DMDc(p=', num2str(p), ')']}, 'Interpreter', 'latex');

    set(gca, 'FontSize', 10, 'LineWidth', 1);
    print('-depsc2', '-r300', fullfile('test', [regexprep(title_label, '[^a-zA-Z0-9]', ''), '.eps']));
end
