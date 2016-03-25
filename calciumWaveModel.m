% This script solves the PDE system for the model proposed in my candidacy
% using the solvepde function

function calciumWaveModel

%% Parameter declaration
% Will use a function later to set default parameters
v_PLC         = 1e-2;   	% [uM/s]    Rate of PLC-gamma
gammaShape    = 1;        % [  ]      Shape constant for noise term
K_Ca       	= 0.2;    	% [uM]      Half-saturation constant for calcium activation of IP3R
K_IP3         = 0.3;      % [uM]      Half-saturation constant for IP3 activation of IP3R
K_i           = 0.09;    	% [uM]      Half-saturation constant for calcium inhibition of IP3R
k_1           = 0.0004;	% [1/s]     Rate constant of calcium leak from ER
k_2           = 0.08;     % [1/s]     Rate constant of calcium release through IP3
k_SERCA         = 0.02;     % [uM]      Half-saturation constant for SERCA pump
v_SERCA           = 0.08;     % [uM/s]    Maximal rate for SERCA pump
k_i           = 13.3;   	% [1/s]     Rate constant of IP3R inactivation
k_deg           = 0.08;     % [1/s]     Rate constant of IP3 degradation     
k_EC      = 0.035;    % [uM]      Saturation constant for calcium transfer with media
v_EC    = 0.51;     % [1/s]     Rate constant for calcium transfer with media

% Physical parameters
D_IP3         = 280;      % [uM/s^2]  Diffusion coefficient of IP3	
D_Ca          = 20;       % [uM/s^2]	Effective diffusion coefficient of calcium
P_IP3         = 1;     	% [uM/s]    Gap-junctional permeability of IP3
P_Ca          = 0.01;     % [uM/s]    Effective gap-junctional permeability of calcium
beta          = 20;       %           Ratio of the effective volumes for calcium of cytoplasm and ER

% Intial conditions
Ca_0	= 0.05;    	% [uM]      Initial concentration of Ca2+ in cytoplasm
CaER_0	= 80;     	% [uM]      Initial concentration of Ca2+ in ER
IP3_0 	= 0.15;    	% [uM]      Initial concentration of IP3
IP3R_0 	= 1;      	% [uM]      Initial active concentration of IP3R

% Domain parameters
totalTime         = 1e4;  % [s]   Total simulation time
cellSize          = 10;	% [um]  Size of cell
cellNumber        = 30;	% [#]	Number of cells in sheet
cellRows          = 30;	% [#]	Number of cells in x direction of sheet
dx = 1;
outputFrames = 100;

%% Initialization
% Homogenize diffusivity
Deff_Ca = 1./(1./D_Ca + 1./(P_Ca .* cellSize));
Deff_IP3 = 1./(1./D_IP3 + 1./(P_IP3 .* cellSize));

% Calculate domain size
L_x = cellSize*cellNumber;
L_y = cellSize*cellRows;

% Lattice points
ts = linspace(0,totalTime,outputFrames);

%% Setup model geometry
% Create empty pde model
model = createpde();

% Decompose a square domain into minimal regions
dl = decsg([3;4;0;L_x;L_x;0;0;0;L_y;L_y],'SQ1',[83;81;49]); 

% Add the decomposed geometry to the PDE model
geometryFromEdges(model, dl);

% Apply Neumann boundary conditions
applyBoundaryCondition(model,'edge',1:model.Geometry.NumEdges,'g',0,'q',0);

% Apply initial conditions
setInitialConditions(model, [Ca_0; IP3_0; CaER_0; IP3R_0], [0;0;0;0]);

%% Describe system of PDE's in m, d, c, a, f form
% No second derivatives with respect to t
m = 0;

% Describes each equation as a generation rate (du/dt = ...)
d = eye(4);
 
% Adds diffusion terms (del2(u))
c = [Deff_Ca ,0     ,0,0; ...
     0    ,Deff_IP3 ,0,0; ...
     0    ,0     ,0,0; ...
     0    ,0     ,0,0];
 
% Adds decay terms (-ku)
a = [k_EC ,0     ,0,0; ...
     0    ,k_deg ,0,0; ...
     0    ,0     ,0,0; ...
     0    ,0     ,0,k_i];

%% Solve system
specifyCoefficients(model,'m',m,'d',d,'c',c,'a',a,'f',@pdeFluxEquations);
generateMesh(model);
results = solvepde(model, ts);
u = results.NodalSolution;
pdeplot3D(model,'colormapdata',u)

function f = pdeFluxEquations(~, state)
% Calculate squares of C and I to reduce computation
C2 = state.u(1,:).^2;
I2 = state.u(2,:).^2;

% Calculate flux terms to reduce compuatation
J_flux = (k_1+k_2.*state.u(4,:).*C2.*I2) ...
    .*(state.u(3,:)-state.u(1,:)) ...
    ./(K_Ca^2+C2) ...
    ./(K_IP3^2+I2);
J_SERCA = v_SERCA .* C2 ./ (K_SERCA^2 + C2);

% Solve for remaining terms
f(1,:) = J_flux - J_SERCA + v_EC;
f(2,:) = v_gen;
f(3,:) = -beta*(J_flux - J_SERCA);
f(4,:) = k_i*(K_i^2./(K_i^2+C2));


