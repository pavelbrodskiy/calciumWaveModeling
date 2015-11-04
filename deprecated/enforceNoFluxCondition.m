function [ A ] = enforceNoFluxCondition( A )
    A(1) = A(2);
    A(end) = A((end-1));
end

