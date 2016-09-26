% This script solves the PDE system for the model proposed in my candidacy
% using the solvepde function
% Pavel A Brodskiy 
% 03/26/16
% Zartman Lab 
% University of Notre Dame

function calciumWaveModel2
% ALL PARAMETERS SHOULD BE IN PRODUCTS OF THE FOLLOWING UNITS:
% time:             [s]
% length:           [um]
% concentration:    [uM]
% particle number:  [1e-18 moles] = [602214 particles]

close all
tic

%% Parameter declaration
% Will use a function later to set default parameters
% Biochemical Parameters
v_PLC   	= 0e-4;   	%1e-2;[uM/s]Mean Rate of PLC-gamma
K_Ca       	= 0.2;    	%0.2 [uM]	Half-saturation constant for calcium activation of IP3R
K_IP3     	= 0.3;      % [uM]      Half-saturation constant for IP3 activation of IP3R
K_i         = 0.09;    	% [uM]      Half-saturation constant for calcium inhibition of IP3R
k_1      	= 4e-4;     % [1/s]     Rate constant of calcium leak from ER
k_2     	= 0.08;     % [1/s]     Rate constant of calcium release through IP3
k_SERCA   	= 0.02;     % [uM]      Half-saturation constant for SERCA pump
v_SERCA    	= 0.08;     %0.08 [uM/s	Maximal rate for SERCA pump
k_i       	= 13.3;   	% [1/s]     Rate constant of IP3R inactivation
k_deg      	= 0.08;     % [1/s]     Rate constant of IP3 degradation     
k_EC        = 1e-3;    %0.035 [uM]     Saturation constant for calcium transfer with media
v_EC        = 1e-3;     % [1/s]     Rate constant for calcium transfer with media

% Physical parameters
D_IP3     	= 280;      % [uM/s^2]  Diffusion coefficient of IP3	
D_Ca       	= 20;       % [uM/s^2]	Effective diffusion coefficient of calcium
P_IP3     	= 1;     	% [uM/s]    Gap-junctional permeability of IP3
P_Ca      	= 0.01;     % [uM/s]    Effective gap-junctional permeability of calcium
beta      	= 20;       %           Ratio of the effective volumes for calcium of cytoplasm and ER

% Intial conditions
Ca_0        = 0.05;    	% [uM]      Initial concentration of Ca2+ in cytoplasm
p.CaER_0	= 80;     	% [uM]      Initial concentration of Ca2+ in ER
IP3_0       = 0;    	% [uM]      Initial concentration of IP3
IP3R_0      = 1;      	% [uM]      Initial active concentration of IP3R

% Domain parameters
totalTime 	= 200;      % [s]       Total simulation time
cellSize  	= 10;       % [um]      Size of cell
cellNumber	= 30;       % [#]       Number of cells in sheet
cellRows  	= 30;       % [#]       Number of cells in x direction of sheet
outputFrames= 100;
noiseParm   = 50;

% Parameters to be passed to f function
p.k_1       = k_1;
p.k_2       = k_2;
p.K_Ca      = K_Ca;
p.K_IP3     = K_IP3;
p.v_SERCA   = v_SERCA;
p.k_SERCA   = k_SERCA;
p.beta      = beta;
p.k_i       = k_i;
p.K_i       = K_i;
p.v_EC      = v_EC;
p.k_EC      = k_EC;
% Random seeds for v_PLC
% p.a_1       = rand(1,noiseParm*2)*1e-1; % Period
% p.a_2       = rand(1,noiseParm*2) * pi .* p.a_1; % Lag
% p.a_3       = rand(1,noiseParm*2)*0+1; % Amplitude
% p.v_PLC     = v_PLC;
% p.noiseParm = noiseParm;
% for k = 1:16
% % for k = 1:16
p.a_1       = rand(1,noiseParm*2)*0.2; % Period
p.a_2       = rand(1,noiseParm*2) * 2* pi; %.* p.a_1; % Lag
p.a_3       = rand(1,noiseParm*2)*0+1; % Amplitude

p.noiseParm = noiseParm;
p.offset = 40;
%     subplot(4,4,k)
% rateMatrix=stochasticMatrix(1:1:300, 1:1:300, p);
p.stoch = @(region)stochasticMatrix(region.x, region.y, p);

% Calculate domain size
L_x = cellSize*cellNumber;
L_y = cellSize*cellRows;

% Normalize PLC generation to mean

stocDist = stochasticMatrix(rand(1,1e5)*L_x, rand(1,1e5)*L_y, p);
% % subplot(4,4,k)
% % hist(stocDist,100)
p.v_PLC     = v_PLC / mean(stocDist); %,std(stocDist)]
% % end
% end
%% Initialization
% Homogenize diffusivity
Deff_Ca = 1./(1./D_Ca + 1./(P_Ca .* cellSize));
Deff_IP3 = 1./(1./D_IP3 + 1./(P_IP3 .* cellSize));


% Lattice points
ts = linspace(0,totalTime,outputFrames);

%% Setup model geometry
% Create empty pde model
model = createpde(3);

% Decompose a square domain into minimal regions
dl = decsg([3;4;0;L_x;L_x;0;0;0;L_y;L_y],'SQ1',[83;81;49]); 

% Add the decomposed geometry to the PDE model
geometryFromEdges(model, dl);

% Apply Neumann boundary conditions
applyBoundaryCondition(model,'edge',1:model.Geometry.NumEdges,'g',zeros(1,3),'q',zeros(1,9));

%% Describe system of PDE's in m, d, c, a, f form
% No second derivatives with respect to t
m = 0;

% Describes each equation as a generation rate (du/dt = ...)
d = [1, 1, 1]';
 
% Adds diffusion terms (del2(u))
c = [Deff_Ca ,Deff_IP3, 0]';
 
% Adds decay terms (-ku)
a = [0, k_deg, 0]';
    
f = @(region, state)pdeFluxEquations(region, state, p);

%% Solve system
specifyCoefficients(model,'m',m,'d',d,'c',c,'a',a,'f',f);
setInitialConditions(model, [Ca_0; IP3_0; IP3R_0]);
generateMesh(model);
results = solvepde(model, ts);
u = results.NodalSolution;

%% Plot final solution (2D)
figure(1)
pdeplot(model,'xydata',u)

%% Plot time-course of small RoI
% [Ca_0; IP3_0; CaER_0; IP3R_0]
h2=figure(2)
subplot(2,2,1)
plot(ts, squeeze(u(100,1,:))); 
xlabel('Time (s)')
ylabel('[Ca^2^+] (uM)')
axis([0,max(ts),0,max(u(100,1,:))*1.2])
subplot(2,2,2)
plot(ts, squeeze(u(100,2,:))); 
xlabel('Time (s)')
ylabel('[IP_3] (uM)')
axis([0,max(ts),0,max(u(100,2,:))*1.2+1e-10])
% subplot(2,2,3)
% plot(ts, squeeze(u(100,3,:))); 
% xlabel('Time (s)')
% ylabel('[Ca^2^+_E_R] (mM)')
% axis([0,max(ts),0,max(u(100,3,:))*1.2])
subplot(2,2,4)
plot(ts, squeeze(u(100,3,:))); 
xlabel('Time (s)')
ylabel('IP_3R Activation (a.u.)')
axis([0,max(ts),0,max(u(100,3,:))*1.2])

%% Save log of simulation
saveName = ['simulationLog' strrep(datestr(now()),':','-')];
save([saveName '.mat'])
savefig(h2,[saveName '.fig'],'compact');

function f1 = pdeFluxEquations(region, state, p)
% Calculate squares of C and I to reduce computation
C2 = state.u(1,:).^2;
I2 = state.u(2,:).^2;
%[state.time, toc]

% Calculate flux terms to reduce compuatation
J_flux = (p.k_1+p.k_2.*state.u(3,:).*C2.*I2) ...
    .*(p.CaER_0-state.u(1,:)) ...
    ./(p.K_Ca^2+C2) ...
    ./(p.K_IP3^2+I2);
J_SERCA = p.v_SERCA .* C2 ./ (p.k_SERCA^2 + C2);
J_EC = p.v_EC - p.k_EC * state.u(1,:); 
% Interpolate sin for continuous noise
[sum(J_flux) sum(J_SERCA) sum(J_EC) state.time, toc]
% Solve for remaining terms
f1(1,:) = (J_flux - J_SERCA + J_EC);
f1(2,:) = p.v_PLC * p.stoch(region);
f1(3,:) = -p.k_i*(p.K_i^2./(p.K_i^2+C2)+state.u(3,:)-1);



function stochMat = stochasticMatrix(x, y, p)
a1s = p.a_1(1:p.noiseParm)' * x;
a2s = repmat(p.a_2(1:p.noiseParm), [length(x),1]);
a3s = p.a_1((p.noiseParm+1):(p.noiseParm*2))' * y;
a4s = repmat(p.a_2((p.noiseParm+1):(p.noiseParm*2)), [length(y),1]);
stochMat = sum(sin(a1s + a2s'),1) +... 
    sum(sin(a3s + a4s'),1);
stochMat = stochMat + p.offset;

% a1s = p.a_1(1:p.noiseParm)' * x;
% a2s = repmat(p.a_2(1:p.noiseParm), [length(x),1]);
% a3s = p.a_1((p.noiseParm+1):(p.noiseParm*2))' * y;
% a4s = repmat(p.a_2((p.noiseParm+1):(p.noiseParm*2)), [length(y),1]);
% stochMat = sum(repmat(p.a_3(1:p.noiseParm),[length(y),1])'.*sin(a1s + a2s'),1) +... 
%     sum(repmat(p.a_3((p.noiseParm+1):(p.noiseParm*2)),[length(y),1])'.*sin(a3s + a4s'),1);
% stochMat = stochMat + p.offset;

% a1s = p.a_1(1:p.noiseParm)' * x;
% a2s = repmat(p.a_2(1:p.noiseParm), [length(x),1]);
% a3s = p.a_1((p.noiseParm+1):(p.noiseParm*2))' * y;
% a4s = repmat(p.a_2((p.noiseParm+1):(p.noiseParm*2)), [length(y),1]);
% stochMat = repmat(sum(repmat(p.a_3(1:p.noiseParm),[length(y),1])'.*sin(a1s + a2s'),1), [length(x), 1]) + repmat(sum(repmat(p.a_3((p.noiseParm+1):(p.noiseParm*2)),[length(y),1])'.*sin(a3s + a4s'),1), [length(x), 1])';
% stochMat = stochMat + p.offset;




