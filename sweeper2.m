function [ values, outputFlag, frequency, amplitude, width ] = sweeper2( meanVal, parm, deviation, runs, timeAtStart )
% meanVal       Double      Mean parameter value
% parm          String      Parameter name
% deviation     Double      Percent deviation from mean (min and max)
% runs          Integer     Total number of parameter values to test
% timeAtStart   Time        OPTIONAL: Time of simulation start

if nargin < 5
   timeAtStart = now();
end

if deviation > 1 % Parse percentage
    deviation = deviation / 100;
end


minVal = meanVal * (1 - deviation);
maxVal = meanVal * (1 + deviation);
stepVal = ((maxVal - minVal) / runs);

i = 1;
values = (minVal:stepVal:maxVal);

if deviation == 0
    values = ones(1,runs) * meanVal;
end


outputDirectory = [ 'Simulation Output' filesep strrep(datestr(timeAtStart),':','-') filesep parm ' Sweep' ];
mkdir(outputDirectory);
outputDirectory = ['''' outputDirectory ''''];

for value = values
    %[ outputFlag(i), summaryStatistics(i) ] = calciumWaveSimulation(['$P' parm], value, 'outputDirectory', outputDirectory);
    [ ~, ~ ] = calciumWaveSimulation([parm], value, 'outputDirectory', outputDirectory);
    i = i + 1;
    disp(['Completed Run ' num2str(i) ' of ' num2str(length(values))]);
end

end

