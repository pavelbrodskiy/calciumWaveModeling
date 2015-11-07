function [ p ] = defaultSettings()
p.outputDirectory   = 'Simulation Output';
p.outputModes       = [3];  % 1 - Analysis, 2 - Final Image, 3 - Realtime Plots, 4 - Video
p.totalTime         = 5e3;  % [s]   Total simulation time
p.cellSize          = 10;	% [um]  Size of cell
p.cellNumber        = 100;	% [#]	Number of cells in sheet
p.cellRows          = 1;	% [#]	Number of cells in x direction of sheet
p.dt                = 5e-2;	% [s]   Simulation timestep
p.dx                = 1;	% [um]  Simulation spatial stepsize
p.outFrames         = 1e3;	% [#]	Number of output frames
p.maxIP3            = 5;	% [uM]  Maximum allowable IP3
p.outputStart       = 1;    % [#]	First frame to be outputted

p.CaBound           = [0, 1];
p.CaERBound       	= [0, 80];
p.IP3Bound      	= [0,10];
p.IP3RBound         = [0, 1];
p.boundCondition    = 'per';% [str] Boundary conditions ('per', 'noflux')
end

