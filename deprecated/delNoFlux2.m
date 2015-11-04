function [ A ] = delNoFlux2( A, dx )

	left  = cat(1, A(2:end   ,      :), A(end       ,         :));
    right = cat(1, A(1       ,      :), A(1:(end-1) ,         :));
    down  = cat(2, A(:       ,      1), A(:         , 1:(end-1)));
    up    = cat(2, A(:       ,  2:end), A(:         ,       end));

    A = (left + right + up + down - 4 * A) / dx ^ 2;
    
end

