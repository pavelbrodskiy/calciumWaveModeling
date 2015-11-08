function [ p ] = defaultSettings()
outputTime = 5e3;

p.outputDirectory   = 'Simulation Output';
p.outputModes       = [2,3];  % 1 - Analysis, 2 - Final Image, 3 - Realtime Plots, 4 - Video
p.totalTime         = 5e4;  % [s]   Total simulation time
p.cellSize          = 4;	% [um]  Size of cell
p.cellNumber        = 30;	% [#]	Number of cells in sheet
p.cellRows          = 1;	% [#]	Number of cells in x direction of sheet
p.dt                = 1e-2;	% [s]   Simulation timestep
p.dx                = (p.dt/5e-2)^2;	% [um]  Simulation spatial stepsize
p.maxIP3            = 1e10;	% [uM]  Maximum allowable IP3
p.kymDim            = [1e3, 1e3];

p.outFrames         = p.totalTime/5;	% [#]	Number of output frames
outputTime          = min([outputTime, p.totalTime]);
p.firstFrame        = round((p.totalTime - outputTime)/dt);


p.outputStart       = 1;    % [#]	First frame to be outputted

p.CaBound           = [0, 0.2];
p.CaERBound       	= [0, 2];
p.IP3Bound      	= [0.1, 0.2];
p.IP3RBound         = [0, 1];
p.boundCondition    = 'noflux';% [str] Boundary conditions ('per', 'noflux')
end

