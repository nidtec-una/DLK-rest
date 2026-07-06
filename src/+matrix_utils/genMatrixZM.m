function [A, b] = genMatrixZM(delta)
    A = [0.5 delta 0 0; 0 1 delta 0; 0 0 1.5 delta; 0 0 0 2];
    b = [0.5 + delta; 1 + delta; 1.5 + delta; 2];
end
