% Implementation of Hofer 2002 astrocyte model
% Pavel Brodskiy, University of Notre Dame, Zartman Lab, 10.21.15
% http://www.jneurosci.org/content/22/12/4850.full

function [ frequency, amplitude, width ] = fHofer2002 (v_7, v_8, P_IP3)
% INPUTS: Input each parameter which is going to be changed as two inputs: 
% first a string identifying the parameter, and then a value.
%
% EXAMPLE: fHofer2002


%% PARAMETER DECLARATION
% Hofer Model Parameters
k_1     = 0.0004;	% [1/s]     Rate constant of calcium leak from ER
k_2     = 0.08;     % [1/s]     Rate constant of calcium release through IP3
v_40    = 0.025;    % [uM/s]    Rate of calcium leak across the plasma membrane
v_41    = 0.2;      % [uM/s]    Maximal rate of activation-dependent calcium influx
k_5     = 0.5;      % [1/s]     Rate constant of calcium extrusion
k_6     = 4;        % [1/s]     Rate constant of IP3R inactivation
k_9     = 0.08;     % [1/s]     Rate constant of IP3 degradation     
K_IP3	= 0.3;      % [uM]      Half-saturation constant for IP3 activation of IP3R
K_a     = 0.2;      % [uM]      Half-saturation constant for calcium activation of IP3R
K_i     = 0.2;      % [uM]      Half-saturation constant for calcium inhibition of IP3R
K_Ca    = 0.3;      % [uM]      Half-saturation constant for calcium activation of PLCd
K_r     = 1;        % [uM]      Half-saturation constant for agonist-dependent calcium entry
D_IP3   = 280;      % [uM/s^2]  Diffusion coefficient of IP3	
D_Ca    = 20;       % [uM/s^2]	Effective diffusion coefficient of calcium
beta    = 20;       %           Ratio of the effective volumes for calcium of cytoplasm and ER

% Tunable Parameters
k_3     = 0.8;          % [1/s]     Rate constant of SERCA pump
% v_7     = 0.08;         % [uM/s]    Maximal rate of PLCd
% v_8     = 1e-2;         % [uM/s]    Rate of PLCb
% P_IP3   = 1;            % [uM/s]    Gap-junctional permeability of IP3
P_Ca    = 0.01*P_IP3;   % [uM/s]    Effective gap-junctional permeability of calcium
alpha_0 = 1.0;          
Kappa_g = 1.0;

%% MODEL PARAMETER DECLARATION
totalTime 	= 2500;     % [s]   Total simulation time
cellSize    = 10;    	% [um]  Size of cell
cellNumber  = 100;    	%       Number of cells in sheet
dt          = 0.05;  	% [s]   Simulation timestep
dx          = 1;        % [um]  Simulation spatial stepsize
outFrames   = 50;    	%       Number of output frames
noiseSigma  = 0.3;      % [1/s] Noise per second

Ca_cyto_0   = 0;        % [uM]  Initial concentration of IP3
Ca_ER_0     = 80;     	% [uM]  Initial concentration of Ca2+ in cytoplasm
IP3_0       = 0;        % [uM]  Initial concentration of Ca2+ in ER
IP3R_0      = 1;        % [uM]  Initial active concentration of IP3R

domainSize = cellNumber * cellSize;
xSize   = 1;
ySize   = domainSize / dx;

IP3_spike = zeros(xSize, ySize);
IP3_spike(uint32(0.4*ySize):uint32(0.6*ySize)) = 3;   	% [uM]  Initial concentration of Ca2+ in ER

stimTime        = 0;
bla = 1;
%% INITIALIZATION

frequency = NaN;
    amplitude = NaN;
    width = NaN;
% Solve for steady state


% Convert to modeling units if needed (uM, um, s, 1e-18 moles)
% Calculate squares of constants

% Initialize variables


CaC     = ones(xSize, ySize) * Ca_cyto_0;
CaER    = ones(xSize, ySize) * Ca_ER_0;
IP3     = ones(xSize, ySize) * IP3_0;
IP3R    = ones(xSize, ySize) * IP3R_0;

v_PLCb  = v_8;
xs      = dx:dx:domainSize;

totalFrames =uint32((totalTime / dt));
framesPerOutput = uint32(totalFrames / outFrames);
frame   = 0;

outputMax = 0.05;

%% BUILD CELL BOUNDARIES

% for i = 2:cellSize/dx:domainSize/dx
%    boundary(int32(i)) = 1; 
% end
% 
% boundary(uint32(domainSize/dx)) = 0;
% 
% boundleft = logical(boundary);
% boundleftleft = logical(boundleft);
% boundright = logical(boundleft);
% boundleftleft(1) = [];
% boundleftleft = [boundleftleft 0];
% boundright = [0 boundright];
% boundright(end) = [];
% boundrightright = logical([0 boundright]);
% boundrightright(end) = [];
% boundleftleft = logical(boundleftleft);
% boundright = logical(boundright);
% doublebound = double(~(boundright | boundleft));

D_Ca = 1/(1/D_Ca + 1/(P_Ca * cellSize))/bla;
D_IP3 = 1/(1/D_IP3 + 1/(P_IP3 * cellSize))/bla;

% CaC_plot = zeros(1, totalFrames);
% CaER_plot = CaC_plot;
% IP3_plot = CaC_plot;
% IP3R_plot = CaC_plot;

CaC_plot = zeros(ySize, totalFrames);
% CaER_plot = CaC_plot;
% IP3_plot = CaC_plot;
% IP3R_plot = CaC_plot;


%% SIMULATION LOOP

for t = 0:dt:totalTime
    if min([CaC CaER IP3 IP3R]) < 0
        return
        plotConcentrations( xs, CaC, CaER, IP3, IP3R, domainSize )
        error('Negative Concentration');      
    end
    
    % Spike IP3
    if t < stimTime
        IP3 = IP3 + IP3_spike * dt;
    end
    %     close all
%     plot(IP3)
    
    % Introduce random noise
    IP3     = IP3 .* normrnd(1,noiseSigma*dt,[xSize, ySize]);
    IP3(IP3<0) = 0;
    
    % Calculate squares for speed
    C2      = CaC  .^ 2;
    I2      = IP3  .^ 2;
    R2      = IP3R .^ 2;
    
    % Calculate rates
    v_rel   = (k_1 + k_2*IP3R.*C2.*I2./(K_a^2+C2)./(K_IP3^2+I2)).*(CaER-CaC);
    v_PLCd  = v_7*C2./(K_Ca^2+C2);
    v_SERCA = k_3*CaC;
    v_out   = k_5*CaC;
    v_deg   = k_9*IP3;
    v_in    = v_40+v_41*I2./(K_r^2+I2);
    dIP3Rdt = k_6*(K_i^2./(K_i^2+C2)-IP3R);
    
    % Calculate laplacian of Ca and IP3
    dif_Ca  = D_Ca * del2Center(CaC, dx); %.* doublebound;
    dif_IP3 = D_IP3 * del2Center(IP3, dx); %.* doublebound;
    
%     % Calculate cell permiability flux
%     jmp_Ca  = P_Ca*(gradient(CaC,dx).*double(boundleft)-gradient(CaC,dx).*double(boundright)); 
%     jmp_IP3 = P_IP3*(gradient(CaC,dx).*double(boundleft)-gradient(IP3,dx).*double(boundright));
    
    % Solve for partial from rates and laplacian
    dCaCdt  = v_rel - v_SERCA + v_in - v_out + dif_Ca;
    dCaERdt = beta * (v_SERCA - v_rel);
    dIP3dt  = v_PLCb + v_PLCd - v_deg + dif_IP3;
    
    % Update concentrations with forward euler method
    CaC     = CaC + dCaCdt * dt;
    CaER    = CaER + dCaERdt * dt;
    IP3     = IP3 + dIP3dt * dt;
    IP3R    = IP3R + dIP3Rdt * dt;
    
    IP3 (IP3>5) = 5;
    
%     % Apply no-flux boundary condition
%     IP3(boundleftleft)      = (IP3(boundleft) + IP3(boundleftleft))/2;
%     IP3(boundrightright)    = (IP3(boundright) + IP3(boundrightright))/2;
%     IP3(boundright)         = IP3(boundrightright);
%     IP3(boundleft)          = IP3(boundleftleft);
%     
%     CaC(boundleftleft)      = (CaC(boundleft) + CaC(boundleftleft))/2;
%     CaC(boundrightright)    = (CaC(boundright) + CaC(boundrightright))/2;
%     CaC(boundright)         = CaC(boundrightright);
%     CaC(boundleft)          = CaC(boundleftleft);
    
    % Output the frame
     if mod(frame, framesPerOutput) ==0
%         
%         %plotConcentrations( xs, CaC, CaER, IP3, IP3R, domainSize );
         disp([num2str(uint16(100*t/totalTime)) '%']);
% %         outputMax = max([outputMax mean(IP3(:)) mean(CaC(:))]);
% %         plot(xs, IP3, xs, CaC);
% %         axis([dx, domainSize, 0, outputMax]);
%         drawnow
     end
    
    frame = frame + 1;
    
%     CaC_plot(frame) = mean(CaC);
%     CaER_plot(frame) = mean(CaER);
%     IP3_plot(frame) = mean(IP3);
%     IP3R_plot(frame) = mean(IP3R);
    
    CaC_plot(:,frame) = CaC;
    %CaER_plot(:,frame) = CaER;
    %IP3_plot(:,frame) = IP3;
    %IP3R_plot(:,frame) = IP3R;

    
end

%% OUTPUT

%close all
%showFig(CaC_plot(:,5000:end)',['v8 = ' num2str(v_8)]);
[ frequency, amplitude, width ] = analyzeWaveOutput( CaC_plot(:,5000:end), dx, dt );
