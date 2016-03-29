% Calcium wave model for Drosophila wing imaginal disc, based roughly on
% Hofer 2002 astrocyte model
%
% Pavel Brodskiy, University of Notre Dame, Zartman Lab, 10.27.15
% Hofer Paper: http://www.jneurosci.org/content/22/12/4850.full
% Narciso Paper: http://iopscience.iop.org/article/10.1088/1478-3975/12/5/056005

function [ outputFlag, summaryStatistics ] = calciumWaveSimulation (varargin)
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

addpath('downloadedDependencies','customDependencies');

tic
close all

%% PARAMETER DECLARATION
p1 = defaultParameters(1);
p2 = defaultSettings();

p = catstruct(p1, p2);

%% INPUT PARSING

% Initialize simulation description variables
domainSizey         = p.cellNumber * p.cellSize;
ySize               = round(domainSizey / p.dx);

if p.cellRows > 1
    domainSizex     = p.cellRows * p.cellSize;
    xSize           = round(domainSizex / p.dx);
    dimensions      = 2;
else
    xSize           = 1;
    dimensions      = 1;
end

[ outputFlag, fileName, p ] = arginParser( p, [xSize, ySize], varargin );

%% INITIALIZATION

% In future version:
% Solve for steady state, convert to modeling units if needed (uM, um, s,
% 1e-18 moles), calculate squares of constants, initialize variables


totalFrames         = uint32((p.totalTime / p.dt));
framesPerOutput     = uint32(totalFrames / p.outFrames);
frame               = 0;
xs                  = p.dx:p.dx:domainSizey;

% Initialize simulation variables
CaC                 = ones(xSize, ySize) * p.Ca_cyto_0;
CaER                = ones(xSize, ySize) * p.Ca_ER_0;
IP3                 = ones(xSize, ySize) * p.IP3_0;
IP3R                = ones(xSize, ySize) * p.IP3R_0;

% Initialize output variables
summaryStatistics.blankField = NaN;
if any(p.outputModes==4)
    aviobj = VideoWriter(['video ' fileName '.avi']);
    open(aviobj);
end

if any(p.outputModes==1)||any(p.outputModes==2)||any(p.outputModes==4)
    if dimensions == 1
        CaC_plot = zeros(ySize, p.kymDim);
    elseif dimensions == 2
        CaC_plot2D = zeros(xSize, ySize, p.outFrames);
    else
        disp('Simulation dimensions unknown');
        outputFlag = 6;
        return
    end
end

% Initialize random pulse variables
timeUntilPulse = exprnd(p.pulseTimeCon(0));
IP3PulseCoords = round(p.IP3Extent / p.dx + 0.5);

%% HOMOGENIZE DIFFUSIVITY

Deff_Ca = 1./(1./p.D_Ca + 1./(p.P_Ca .* p.cellSize));
Deff_IP3 = 1./(1./p.D_IP3 + 1./(p.P_IP3 .* p.cellSize));

%randomNumbers = zeros(xSize, ySize);
plotIterator = 1;
outputTime = p.firstTime;    

%meanProduction = p.v_PLC;
theta = p.v_PLC ./ p.gammaShape;
%v_PLC = p.v_PLC;

%% SIMULATION LOOP
    
for t = 0:p.dt:p.totalTime
    
    % Implement random flashes
    if timeUntilPulse <= 0
        [xPulse, yPulse] = pulseCoordinates(xSize, ySize, IP3PulseCoords, p);
        IP3(xPulse, yPulse) = IP3(xPulse, yPulse) + p.IP3Pulse;
        timeUntilPulse = exprnd(p.pulseTimeConstant(t));
        disp('Pulse Occured');
    else
        timeUntilPulse = timeUntilPulse - p.dt;
    end
    
    % Bookeeping
    outputFlag = evaluateFlags( CaC, CaER, IP3, IP3R, p.maxIP3 ) > 0;
    if outputFlag > 0 
        return 
    end
    
    % Calculate squares for speed
    C2          = CaC .^ 2;
    I2          = IP3 .^ 2;
    
    % Calculate rates
    v_rel   = (p.k_1 + p.k_2.*IP3R.*C2.*I2./(p.K_Ca.^2+C2)./(p.K_IP3.^2+I2)).*(CaER-CaC);
    v_SERCA = (p.gam .* C2) ./ (p.k_gam.^2 + C2);
    v_deg   = p.k_9.*IP3;
    v_media    =  p.P_Ca_media .* (p.Ca_media - CaC);
    dIP3Rdt = p.k_6.*(p.K_i.^2./(p.K_i.^2+C2)-IP3R);
    v_PLC = gamrnd(p.gammaShape,theta, [xSize, ySize]);
    
    % Calculate laplacian of Ca and IP3 with no-flux boundary conditions
    switch p.boundCondition
        case 'per'
            if dimensions == 1
                dif_Ca	= Deff_Ca .* del2Periodic1D(CaC, p.dx);
                dif_IP3	= Deff_IP3 .* del2Periodic1D(IP3, p.dx);
            elseif dimensions == 2
                dif_Ca	= Deff_Ca .* del2Periodic2D(CaC, p.dx);
                dif_IP3	= Deff_IP3 .* del2Periodic2D(IP3, p.dx);
            end
        case 'noflux'
            if dimensions == 1
                dif_Ca	= Deff_Ca .* del2NoFlux1D(CaC, p.dx);
                dif_IP3	= Deff_IP3 .* del2NoFlux1D(IP3, p.dx);
            elseif dimensions == 2
                dif_Ca	= Deff_Ca .* del2NoFlux2D(CaC, p.dx);
                dif_IP3	= Deff_IP3 .* del2NoFlux2D(IP3, p.dx);
            end
        otherwise
            error('incorrect boundary')
    end
    
    % Solve for partial from rates and laplacian
    dCaCdt  = v_rel - v_SERCA + v_media + dif_Ca;
    dCaERdt = p.beta .* (v_SERCA - v_rel);
    dIP3dt  = v_PLC - v_deg + dif_IP3;
    
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
                %plotOneConcentration( xs, CaC, domainSizey, p );
                plotConcentrations( xs, CaC, CaER, IP3, IP3R, domainSizey, p )
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
            CaC_plot2D(:,:,frame) = CaC;
        end
    end
    
    % Record concentration of calcium
    frame = frame + 1;
    if t>=outputTime&&(any(p.outputModes==1)||any(p.outputModes==2)||any(p.outputModes==4))
        if dimensions == 1
            CaC_plot(:,plotIterator) = CaC;
            plotIterator = plotIterator + 1;
            outputTime = t + p.outputInterval;
        end
    end
end

%% OUTPUT
% Use the concentration record over time to automatically generate
% simulation output

close all

if any(p.outputModes==1) % Calculate summary statistics
    disp('Calculating simulation properties');
    [ summaryStatistics.frequency, summaryStatistics.amplitude, summaryStatistics.width ] = analyzeWaveOutput( CaC_plot(:,p.outputStart:end), p.dx, p.dt );
end
if any(p.outputModes==2) % Generate kymograph
    disp(['Saving kymograph: ' fileName]);
    showFig(CaC_plot',[p.outputDirectory filesep fileName],p);
end
if any(p.outputModes==4) % I broke this, this makes a video over time
    error('THIS HAS NOT BEEN WRITTEN');
    writeVideo(aviobj, getframe(gcf))
    disp(['Saving video: ' fileName]);
    close(aviobj)
end

toc
summaryStatistics.runTime = toc;

