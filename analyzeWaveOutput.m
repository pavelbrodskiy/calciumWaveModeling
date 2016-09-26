function [ frequency, amplitude, width ] = analyzeWaveOutput( twoDimensionalSignal, dx, dt )
% Input an M x N matrix where M is the x dimension and N is the time
% dimension. The waves are decomposed into individual x positions, then
% amplitude and frequency is outputted as a function of postion.
%
% ALGORITHM
% 1. Decompose signal into bins in space and time. Could use gaussian blur
%       in space to reduce effect of noise.
% 2. Obtain the time at which peaks occur for each spatial coordinate. Can
%       use FFT to obtain frequency or peaks function.
% 3. Use peaks function with minimum time between peaks as equal to 2/3 of
%       the measured time. Some smoothing may be needed to prevent false
%       positives, or some minimum peak height or width.
% 4. Frequency is the median of 1/(time between each wave). Amplitude is
%       the median of (peak value) - (average minimum value). Wave width
%       is the median of width (from peaks output). Wave speed is the
%       median of the slope of front position and time. Regime outputs
%       what the "qualitative" behavior of the simulation is as a function 
%       of time. 
%
% Frequency     [double]    	1/(time between each wave)
% Amplitude     [double]        (peak value) - (average minimum value)
% Wave width  	[double]
% Wave speed 	[double]        
% Regime        [int array 1xT]	1 = waves, 2 = pulses, 3 = global 
%                                   oscilation, 4 = steady state over time

[totalX, totalT] = size(twoDimensionalSignal);

for i = totalX:-1:1
   signal = twoDimensionalSignal(i,:);
   [pks,locs,w,p] = findpeaks(signal,'MinPeakDistance',2000);
   amplitude(i) = mean(pks);
   width = mean(w);
   
   distances = [locs 0] - [0 locs];
   frequency(i) = 1/(mean(distances(2:(end-2)))*dt);
end

amplitude = median(amplitude);
frequency = median(frequency);
width = median(width);
end

