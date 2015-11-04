% Calcium wave model for Drosophila wing imaginal disc, based roughly on
% Hofer 2002 astrocyte model
%
% Pavel Brodskiy, University of Notre Dame, Zartman Lab, 10.27.15
% Hofer Paper: http://www.jneurosci.org/content/22/12/4850.full
% Narciso Paper: http://iopscience.iop.org/article/10.1088/1478-3975/12/5/056005

function [ outputFlag, frequency, amplitude, width ] = calciumWaveSimulation4_D (varargin)
% INPUTS: Input each parameter which is going to be changed as two inputs:
% first a string identifying the parameter, and then a value. Use $A or $P
% to selectivly change parameters in the A or P compartment respectively.
%
% EXAMPLE: [ outputFlag, frequency, amplitude, width ] = calciumWaveSimulation2 ('v_7', 0.02);
% EXAMPLE: [ outputFlag, frequency, amplitude, width ] = calciumWaveSimulation2 ('$Pv_7', 0.02);
%
% OUTPUTS: The simulation will run once, then output a frequency, amputude,
% and a width at half height for each position in the domain. The
% function will also save a chymograph of the simulation with a name
% derived from the modified inputs.
% 
% Frequency (1/period)
% Amplitude
% Wave width
% Wave speed


tic
close all

%% PARAMETER DECLARATION
p1 = defaultParameters();
p2 = defaultSettings();

p = catstruct(p1, p2);

%% INPUT PARSING

[ outputFlag, fileName, p ] = arginParser( p, varargin );


%% INITIALIZATION

% In future version:
% Solve for steady state, convert to modeling units if needed (uM, um, s,
% 1e-18 moles), calculate squares of constants, initialize variables

domainSizey         = p.cellNumber * p.cellSize;
ySize               = domainSizey / p.dx;
if p.cellRows > 1
    domainSizex     = p.cellRows * p.cellSize;
    xSize           = domainSizex / p.dx;
    dimensions      = 2;
else
    xSize           = 1;
    dimensions      = 1;
end

frequency   = NaN;
amplitude   = NaN;
width       = NaN;

CaC         = ones(xSize, ySize) * p.Ca_cyto_0;
CaER        = ones(xSize, ySize) * p.Ca_ER_0;
IP3         = ones(xSize, ySize) * p.IP3_0;
IP3R        = ones(xSize, ySize) * p.IP3R_0;

xs          = p.dx:p.dx:domainSizey;

totalFrames = uint32((p.totalTime / p.dt));
framesPerOutput = uint32(totalFrames / p.outFrames);
frame       = 0;

if any(p.outputModes==4)
    aviobj = VideoWriter(['video ' fileName '.avi']);
    open(aviobj);
end

addpath('downloadedDependencies','customDependencies');

%% HOMOGENIZE DIFFUSIVITY

Deff_Ca = 1./(1./p.D_Ca + 1./(p.P_Ca .* p.cellSize));
Deff_IP3 = 1./(1./p.D_IP3 + 1./(p.P_IP3 .* p.cellSize));

if any(p.outputModes==1)||any(p.outputModes==2)||any(p.outputModes==4)
    if dimensions == 1
        CaC_plot = zeros(xSize, ySize, totalFrames);
    elseif dimensions == 2
        CaC_plot = zeros(xSize, ySize, p.outFrames);
    else
        disp('Simulation dimensions unknown');
        outputFlag = 6;
        return
    end
end

%% SIMULATION LOOP

for t = 0:p.dt:p.totalTime
    
    IP3(IP3 > p.maxIP3) = p.maxIP3;
    
    outputFlag = evaluateFlags( CaC, CaER, IP3, IP3R, p.maxIP3 ) > 0;
    if outputFlag > 0 
        return 
    end
    
    % Introduce random noise
    IP3         = IP3 .* normrnd(1,p.noiseSigma * p.dt, [xSize, ySize]);
    
    % Calculate squares for speed
    C2          = CaC .^ 2;
    I2          = IP3 .^ 2;
    
    % Calculate rates
    v_PLC       =  p.PLC(C2);
    v_rel       = (p.k_1 + p.k_2 .* IP3R .* C2 .* I2 ./ (p.K_a.^2 + C2)./(p.K_IP3.^2 + I2)) .* (CaER - CaC);
    v_SERCA     =  p.k_3 .* CaC;
    v_out       =  p.k_5 .* CaC;
    v_deg       =  p.k_9 .* IP3;
    v_in        =  p.v_40 + p.v_41 .* I2 ./ (p.K_r.^2 + I2);
    dIP3Rdt     =  p.k_6 .* (p.K_i.^2 ./ (p.K_i.^2 + C2) - IP3R);
    
    % Calculate laplacian of Ca and IP3 with no-flux boundary conditions
    if dimensions == 1
        dif_Ca	= Deff_Ca .* del2Center1D(CaC, p.dx);
        dif_IP3	= Deff_IP3 .* del2Center1D(IP3, p.dx);
    elseif dimensions == 2
        dif_Ca	= Deff_Ca .* del2Center2D(CaC, p.dx);
        dif_IP3	= Deff_IP3 .* del2Center2D(IP3, p.dx);
    end
    
    % Solve for partial from rates and laplacian
    dCaCdt      = v_rel - v_SERCA + v_in - v_out + dif_Ca;
    dCaERdt     = p.beta .* (v_SERCA - v_rel);
    dIP3dt      = v_PLC - v_deg + dif_IP3;
    
    % Update concentrations with forward euler method
    CaC         = CaC  + dCaCdt  * p.dt;
    CaER        = CaER + dCaERdt * p.dt;
    IP3         = IP3  + dIP3dt  * p.dt;
    IP3R        = IP3R + dIP3Rdt * p.dt;
    
    % Output the frame realtime
    if mod(frame, framesPerOutput) == 0
        disp([num2str(100*t/p.totalTime) '% done with simulation']);
        if any(p.outputModes==3) % Check if realtime output is requested
            if dimensions == 1
                plotOneConcentration( xs, CaC, domainSizey );
                drawnow
            elseif dimensions == 2
                imshow(CaC, p.CaBound);
                colormap('jet');
                drawnow
            else
                disp('Simulation dimensions unknown');
                outputFlag = 6;
                return
            end
        end
        if dimensions == 2 && frame
            CaC_plot(:,:,frame) = CaC;
        end
    end
    
    % Record concentration of Calcium
    frame = frame + 1;
    if any(p.outputModes==1)||any(p.outputModes==2)||any(p.outputModes==4)
        if dimensions == 1
            CaC_plot(:,:,frame) = CaC;
        end
    end
end

%% OUTPUT
% Use the concentration record over time to automatically generate
% simulation output

close all

if any(p.outputModes==1) % Calculate summary statistics
    disp('Calculating simulation properties');
    [ frequency, amplitude, width ] = analyzeWaveOutput( CaC_plot(:,p.outputStart:end), p.dx, p.dt );
end
if any(p.outputModes==2) % Generate kymograph
    disp(['Saving chymograph: ' fileName]);
    showFig(CaC_plot(:,p.outputStart:end)',[p.outputDirectory filesep p.fileName]);
end
if any(p.outputModes==4) % I broke this, this makes a video over time
    error('THIS HAS NOT BEEN WRITTEN');
    writeVideo(aviobj, getframe(gcf))
    disp(['Saving video: ' fileName]);
    close(aviobj)
end

toc

