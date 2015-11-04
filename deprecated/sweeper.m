function [ values, outputFlag, frequency, amplitude, width ] = sweeper( meanVal, parm, deviation, runs, timeAtStart )
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

outputDirectory = [ strrep(datestr(timeAtStart),':','-') filesep parm ' Sweep' ];
mkdir(outputDirectory);
outputDirectory = ['''' outputDirectory ''''];

for value = values
    [ outputFlag(i), frequency(i), amplitude(i), width(i) ] = calciumWaveSimulation4(['$P' parm], value, 'outputDirectory', outputDirectory, 'outputStart', 20000);
    i = i + 1;
    disp(['Completed Run ' num2str(i) ' of ' num2str(length(values))]);
end

end

