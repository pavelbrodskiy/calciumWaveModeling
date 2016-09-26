function [ p ] = defaultSettings()
outputTime = 2e3;

% Simulation settings
p.outputDirectory   = 'Simulation Output';
p.outputModes       = [3];  % 1 - Analysis, 2 - Final Image, 3 - Realtime Plots, 4 - Video, 5 - Final Figure
p.totalTime         = 1e4;  % [s]   Total simulation time
p.cellSize          = 10;	% [um]  Size of cell
p.cellNumber        = 30;	% [#]	Number of cells in sheet
p.cellRows          = 1;	% [#]	Number of cells in x direction of sheet
p.dt                = 5e-2;	% [s]   Simulation timestep
p.maxIP3            = 1e10;	% [uM]  Maximum allowable IP3
p.kymDim            = 1e3; 
p.dx                = (p.dt/5e-2)^2;	% [um]  Simulation spatial stepsize

% Output settings
p.outFrames         = p.totalTime/5;	% [#]	Number of output frames
outputTime          = min([outputTime, p.totalTime]);
p.firstTime         = max(0,p.totalTime - outputTime);
outputFrames        = min(round(outputTime/p.dt), p.kymDim);
p.outputInterval    = outputTime/outputFrames;
p.outputStart       = 1;    % [#]	First frame to be outputted

p.CaBound           = [0, 0.2];
p.CaERBound       	= [0, 100];
p.IP3Bound      	= [0.1, 0.2];
p.IP3RBound         = [0, 1];

p.boundCondition    = 'noflux';% [str] Boundary conditions ('per', 'noflux')
end