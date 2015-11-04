function [ A ] = del2Periodic1D( A, dx )

    left = [A(2:end) A(1)];
    right = [A(end) A(1:(end-1))];

    A = (left + right - 2 * A) / dx ^ 2;
    
end

