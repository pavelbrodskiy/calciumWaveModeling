% Calcium wave model for Drosophila wing imaginal disc, based roughly on
% Hofer 2002 astrocyte model
%
% Pavel Brodskiy, University of Notre Dame, Zartman Lab, 10.27.15
% Hofer Paper: http://www.jneurosci.org/content/22/12/4850.full

function [ outputFlag, frequency, amplitude, width ] = calciumWaveSimulation4 (varargin)
% INPUTS: Input each parameter which is going to be changed as two inputs:
% first a string identifying the parameter, and then a value. Use $A or $P
% to selectivly change parameters in the A or P compartment respectively.
%
% EXAMPLE: [ frequency, amplitude, width ] = calciumWaveSimulation2 ('v_7', 0.02);
% EXAMPLE: [ frequency, amplitude, width ] = calciumWaveSimulation2 ('$Pv_7', 0.02);
%
% OUTPUTS: The simulation will run once, then output a frequency, amputude,
% and a width at half height for each position in the domain. The
% function will also save a chymograph of the simulation with a name
% derived from the modified inputs.

tic

%% PARAMETER DECLARATION
% Hofer Model Parameters
k_1     = 0.0004;       % [1/s]     Rate constant of calcium leak from ER
k_2     = 0.08;         % [1/s]     Rate constant of calcium release through IP3
v_40    = 0.025;        % [uM/s]    Rate of calcium leak across the plasma membrane
v_41    = 0.2;          % [uM/s]    Maximal rate of activation-dependent calcium influx
k_5     = 0.5;          % [1/s]     Rate constant of calcium extrusion
k_6     = 4;            % [1/s]     Rate constant of IP3R inactivation
k_9     = 0.08;         % [1/s]     Rate constant of IP3 degradation
K_IP3	= 0.3;          % [uM]      Half-saturation constant for IP3 activation of IP3R
K_a     = 0.2;          % [uM]      Half-saturation constant for calcium activation of IP3R
K_i     = 0.2;          % [uM]      Half-saturation constant for calcium inhibition of IP3R
K_r     = 1;            % [uM]      Half-saturation constant for agonist-dependent calcium entry
D_IP3   = 280;          % [uM/s^2]  Diffusion coefficient of IP3
D_Ca    = 20;           % [uM/s^2]	Diffusion coefficient of calcium
beta    = 20;           %           Ratio of the effective volumes for calcium of cytoplasm and ER

% Tunable Parameters
V_CaInd = 8e-3;
V_Ca    = 8e-2;
K_Ca    = 3e-1;
k_3     = 0.8;          % [1/s]     Rate constant of SERCA pump
P_IP3   = 1;            % [uM/s]    Gap-junctional permeability of IP3
P_Ca    = 0.01*P_IP3;   % [uM/s]    Effective gap-junctional permeability of calcium


%v_7     = 0.08;        % [uM/s]    Maximal rate of PLCd
%v_8     = 1e-2;        % [uM/s]    Rate of PLCb
%K_Ca    = 0.3;         % [uM]      Half-saturation constant for calcium activation of PLCd


% Postulated PLC kinetics
PLC = @(Ca2, V_CaDep, K_CaDep, V_CaInd) V_CaInd + V_CaDep.*Ca2./(K_CaDep.^2+Ca2);

%% MODEL PARAMETER DECLARATION
outputDirectory = 'Simulation Output';

outputModes = [2];   % 1 - Analysis, 2 - Final Image, 3 - Realtime Plots, 4 - Video
totalTime 	= 2500;     % [s]   Total simulation time
cellSize    = 10;    	% [um]  Size of cell
cellNumber  = 100;    	%       Number of cells in sheet
dt          = 0.05;  	% [s]   Simulation timestep
dx          = 1;        % [um]  Simulation spatial stepsize
outFrames   = 50;    	%       Number of output frames
noiseSigma  = 0.1;     % [1/s] Noise per second
%noiseSigma  = 0;        % [1/s] Noise per second
maxIP3      = 5;        % [uM]  Maximum allowable IP3
outputStart = 1;

Ca_cyto_0   = 0;        % [uM]  Initial concentration of IP3
Ca_ER_0     = 80;     	% [uM]  Initial concentration of Ca2+ in cytoplasm
IP3_0       = 0;        % [uM]  Initial concentration of Ca2+ in ER
IP3R_0      = 1;        % [uM]  Initial active concentration of IP3R

domainSize = cellNumber * cellSize;
xSize   = 1;
ySize   = domainSize / dx;

% [uM]  Initial concentration of Ca2+ in ER
IP3_spike = zeros(xSize, ySize);
IP3_spike(:,uint32(0.1*ySize):uint32(0.2*ySize)) = 3e-1;
stimTime = 0;

outputFlag = 0;




%% INPUT PARSING
% Note to future-Pavel PLEASE PLEASE switch to input parser instead of eval

if mod(nargin,2) ~= 0
    disp('Inputs should be in the form (''Variable name'', Parameter Value)');
    outputFlag = 4;
    return
end
for i = 1:2:nargin
    temp = varargin{i};
    if ~(temp(1)=='$') % If the input is not compartment-specific
        eval([temp ' = ' num2str(varargin{i+1}) ';']);
    else % If input is compartment-specific, this is indicated with $
        compartment = ones(xSize, ySize);
        
        if (temp(2)=='P') % Check if P or A
            compartment(:, 1:round(ySize/2)) = 0;
        elseif (temp(2)=='A')
            compartment(:, round(ySize/2):end) = 0;
        else
            error('$ flag used without A or P designation');
        end
        
        % Implement the change only in that compartment
        par = temp(3:end);
        if (exist(par,'var'))
            eval([par ' = compartment .* ' num2str(varargin{i+1}) ' + ~compartment .* ' par ';']);
        else
            disp(['Parameter does not exist: ' par]);
            outputFlag = 5;
            return
        end
    end
end

fileName = '';
for i = 1:2:nargin
    if ~(strcmp(varargin{i},'outputDirectory')||strcmp(varargin{i},'outputStart'))
        fileName = [fileName varargin{i} ' = ' num2str(varargin{i+1},'%.6f') ' '];
    end
end


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

%v_PLCb  = v_8;
xs      = dx:dx:domainSize;

totalFrames =uint32((totalTime / dt));
framesPerOutput = uint32(totalFrames / outFrames);
frame   = 0;

if any(outputModes==3)||any(outputModes==4)
    close all;
end
if any(outputModes==4)
    aviobj = VideoWriter(['video ' fileName '.avi']);
    open(aviobj);
end

%% HOMOGENIZE DIFFUSIVITY

D_Ca = 1./(1./D_Ca + 1./(P_Ca .* cellSize));
D_IP3 = 1./(1./D_IP3 + 1./(P_IP3 .* cellSize));

CaC_plot = zeros(ySize, totalFrames);

%% SIMULATION LOOP

for t = 0:dt:totalTime
    
    IP3(IP3 > maxIP3) = maxIP3;
    
    outputFlag = evaluateFlags( CaC, CaER, IP3, IP3R, maxIP3 ) > 0;
    if outputFlag > 0
        return
    end
    
    % Spike IP3
    if t < stimTime
        IP3 = IP3 + IP3_spike * dt;
    end
    
    % Introduce random noise
    IP3 = IP3 .* normrnd(1,noiseSigma*dt,[xSize, ySize]);
    IP3(IP3<0) = 0;
    
    % Calculate squares for speed
    C2      = CaC  .^ 2;
    I2      = IP3  .^ 2;
    
    % Calculate rates
    v_rel   = (k_1 + k_2.*IP3R.*C2.*I2./(K_a.^2+C2)./(K_IP3.^2+I2)).*(CaER-CaC);
    %v_PLCd  = v_7*C2./(K_Ca^2+C2);
    v_PLC   = PLC(C2, V_Ca, K_Ca, V_CaInd);
    v_SERCA = k_3 .* CaC;
    v_out   = k_5 .* CaC;
    v_deg   = k_9 .* IP3;
    v_in    = v_40 + v_41.*I2./(K_r.^2+I2);
    dIP3Rdt = k_6 .* (K_i.^2./(K_i.^2+C2)-IP3R);
    
    % Calculate laplacian of Ca and IP3 (del2Center assumes no flux BC)
    dif_Ca  = D_Ca .* del2Center(CaC, dx);
    dif_IP3 = D_IP3 .* del2Center(IP3, dx);
    
    % Solve for partial from rates and laplacian
    dCaCdt  = v_rel - v_SERCA + v_in - v_out + dif_Ca;
    dCaERdt = beta .* (v_SERCA - v_rel);
    dIP3dt  = v_PLC - v_deg + dif_IP3;
    
    % Update concentrations with forward euler method
    CaC     = CaC + dCaCdt * dt;
    CaER    = CaER + dCaERdt * dt;
    IP3     = IP3 + dIP3dt * dt;
    IP3R    = IP3R + dIP3Rdt * dt;
    
    % Output the frame
    if mod(frame, framesPerOutput) == 0
        if any(outputModes==3)||any(outputModes==4)
            plotOneConcentration( xs, CaC, domainSize );
            drawnow
        end
        if any(outputModes==4)
            writeVideo(aviobj, getframe(gcf))
        end
        disp([num2str(uint16(100*t/totalTime)) '% done with simulation']);
    end
    
    frame = frame + 1;
    
    CaC_plot(:,frame) = CaC;
    
end

%% OUTPUT

close all

if any(outputModes==2)
    disp(['Saving chymograph: ' fileName]);
    showFig(CaC_plot(:,outputStart:end)',[outputDirectory filesep fileName]);
end
if any(outputModes==1)
    disp('Calculating simulation properties');
    [ frequency, amplitude, width ] = analyzeWaveOutput( CaC_plot(:,outputStart:end), dx, dt );
end
if any(outputModes==4)
    disp(['Saving video: ' fileName]);
    close(aviobj)
end

toc

