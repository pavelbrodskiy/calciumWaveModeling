function [ A ] = del2Center( A, dx )

	left  = cat(1, A(2:end   ,      :), A(1         ,         :));
    right = cat(1, A(end     ,      :), A(1:(end-1) ,         :));
    down  = cat(2, A(:       ,    end), A(:         , 1:(end-1)));
    up    = cat(2, A(:       ,  2:end), A(:         ,         1));

    A = (left + right + up + down - 4 * A) / dx ^ 2;
    
end

