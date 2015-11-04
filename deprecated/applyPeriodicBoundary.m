function [ B ] = applyPeriodicBoundary( A )
    B = [A(end) A A(1)];
end

