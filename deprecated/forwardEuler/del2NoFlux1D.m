function [ A ] = del2NoFlux1D( A, dx )

    left = [A(2:end) A(end)];
    right = [A(1) A(1:(end-1))];

    A = (left + right - 2 * A) / dx ^ 2;
    
end

