% This model represents cells as lattice points (assume hexagonal, planar
% lattice) and solves the calcium model as a series of ODEs connected by
% jump boundary conditions.
%
% Solver is forward euler method

close all
clearvars

%% Parameter Declaration  
% Biochemical Parameters
v_PLC   	= 1e-2;   	%1e-2;[uM/s]Rate of PLC-gamma
K_Ca       	= 0.8;    	%0.2 [uM]	Half-saturation constant for calcium activation of IP3R
K_IP3     	= 0.3;      % [uM]      Half-saturation constant for IP3 activation of IP3R
K_i         = 0.09;    	% [uM]      Half-saturation constant for calcium inhibition of IP3R
k_1      	= 4e-4;     % [1/s]     Rate constant of calcium leak from ER
k_2     	= 0.08;     % [1/s]     Rate constant of calcium release through IP3
k_SERCA   	= 0.02;     % [uM]      Half-saturation constant for SERCA pump
v_SERCA    	= 0.8;      %0.08 [uM/s	Maximal rate for SERCA pump
k_i       	= 13.3;   	% [1/s]     Rate constant of IP3R inactivation
k_deg      	= 0.08;     % [1/s]     Rate constant of IP3 degradation     
k_EC        = 7;        %0.035 [uM]	Saturation constant for calcium transfer with media
v_EC        = 0.51;     % [1/s]     Rate constant for calcium transfer with media
gammaShape 	= 1;        % [  ]      Shape constant for noise term

% Physical parameters
D_IP3     	= 280;      % [uM/s^2]  Diffusion coefficient of IP3	
D_Ca       	= 20;       % [uM/s^2]	Effective diffusion coefficient of calcium
P_IP3     	= 1;     	% [uM/s]    Gap-junctional permeability of IP3
P_Ca      	= 0.01;     % [uM/s]    Effective gap-junctional permeability of calcium
beta      	= 20;       %           Ratio of the effective volumes for calcium of cytoplasm and ER

% Intial conditions
Ca_0        = 0.05;    	% [uM]      Initial concentration of Ca2+ in cytoplasm
CaER_0      = 80;     	% [uM]      Initial concentration of Ca2+ in ER
IP3_0       = 0;    	% [uM]      Initial concentration of IP3
IP3R_0      = 1;      	% [uM]      Initial active concentration of IP3R

% Domain parameters
endTime 	= 50;       % [s]       Total simulation time
cellLength 	= 10;       % [um]      Size of cell
outputFrames= 100;
cellsAcross = 6;
dt = 1e-3;


%% Initialize lattice
dx = sqrt(3) * cellLength / 2;
dy = cellLength / 2;
xPositions = 0:dx:dx*(cellsAcross*12+1);
xPositions = repmat(xPositions, [cellsAcross*12+2,1]);
yPositions = 0:dy:dy*(cellsAcross*12+1);
yPositions = repmat([yPositions; yPositions + dx], [cellsAcross*6+1,1])';

xtemp = xPositions(:);
ytemp = yPositions(:);
distances = sqrt((xtemp - cellsAcross * cellLength - 1).^2 + (ytemp - cellsAcross * cellLength - 1).^2);
ytemp(distances > (cellsAcross * cellLength * 1.05)) = [];
xtemp(distances > (cellsAcross * cellLength * 1.05)) = [];

clear distances
adjacencyMatrix = false(length(xtemp));
for i = 1:length(xtemp)
    distances = sqrt((xtemp - xtemp(i)).^2 + (ytemp - ytemp(i)).^2);
    sort(distances,'descend')
    adjacencyMatrix(i,distances>cellLength*0.99) = 0;
    adjacencyMatrix(i,distances<cellLength*0.99) = 1;
    adjacencyMatrix(i,distances==0) = 0;
end

%% Simulation loop
Ca(1:length(xtemp)) = Ca_0;
CaER(1:length(xtemp)) = CaER_0;
IP3(1:length(xtemp)) = IP3_0;
IP3R(1:length(xtemp)) = IP3R_0;

for t = 0:dt:endTime
    C2 = Ca.^2;
    I2 = IP3.^2;
    
    % Calculate rates
    v_rel       = (p.k_1 + p.k_2.*IP3R.*C2.*I2./(p.K_Ca.^2+C2)./(p.K_IP3.^2+I2)).*(CaER-Ca);
    v_SERCA     = (p.gam .* C2) ./ (p.k_gam.^2 + C2);
    v_deg       = p.k_9.*IP3;
    v_media     =  p.P_Ca_media .* (p.Ca_media - Ca);
    dIP3Rdt     = p.k_6.*(p.K_i.^2./(p.K_i.^2+C2)-IP3R);
    %v_PLC       = gamrnd(p.gammaShape,theta, [xSize, ySize]);
    
    % Solve for partial from rates and laplacian
    dCaCdt      = v_rel - v_SERCA + v_media + dif_Ca;
    dCaERdt     = p.beta .* (v_SERCA - v_rel);
    dIP3dt      = v_PLC - v_deg + dif_IP3;
    
    % Update concentrations with forward euler method
    Ca          = Ca  + dCaCdt  * p.dt;
    CaER        = CaER + dCaERdt * p.dt;
    IP3         = IP3  + dIP3dt  * p.dt;
    IP3R        = IP3R + dIP3Rdt * p.dt;
    
    if mod(t,endTime/100)
        scatter(xtemp, ytemp, C2, 3, 'filled')
    end
end

%% Process output
figure(1)
scatter(xtemp(:),ytemp(:),'filled')
axis('square')
% figure(2)
voronoi(xtemp,ytemp)
hold on
cellID=15;
scatter(xtemp(adjacencyMatrix(cellID,:)),ytemp(adjacencyMatrix(cellID,:)),'filled','red')
hold on
scatter(xtemp(cellID),ytemp(cellID),'filled','green')
