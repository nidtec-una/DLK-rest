function [A, b] = genMatrixZMn(delta,n)
    % Inicializar matriz y vector
    A = zeros(n);
    b = zeros(n,1);
    
    % Construir matriz y vector
    for i = 1:n
        A(i,i) = 0.5*i;           % diagonal principal: 0.5, 1, 1.5, ...
        if i < n
            A(i,i+1) = delta;     % superdiagonal: delta
        end
        if i < n
            b(i) = 0.5*i + delta; % b con delta en todas menos la última
        else
            b(i) = 0.5*i;         % último igual que la diagonal
        end
    end
end
