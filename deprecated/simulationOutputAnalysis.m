% Dim 1 is: run number
% Dim 2 is: parameter tested
% Dim 3 is: parameter value in sequence

frequency1 = frequency;
amplitude1 = amplitude;
width1 = width;

frequency1(frequency1==0) = NaN;
amplitude1(amplitude1==0) = NaN;
width1(width1==0) = NaN;

v7 = 0.06:0.01:0.16;
P_IP3 = [0, 0.1, 0.25, 0.5, 1, 2, 3.5, 5, 7.5, 10];
v8 = [0, 1e-4, 8e-4, 2e-3, 1e-2, 5e-2, 0.1, 0.5, 1, 3];

makePlot( width1, v7, P_IP3, v8 )
