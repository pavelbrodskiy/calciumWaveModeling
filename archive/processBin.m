function stats = processBin (videoBin, settings)
% This function takes a video bin and converts it into a 1D signal, then
% extracts summary statistics.

if nargin < 2
    settings = settings();
end

maxProj = max(videoBin,[],3);
minProj = mean(mean(min(videoBin,[],3)));
signalSize = sum(sum(maxProj > 0));
meanIntensity = squeeze(sum(sum(videoBin - minProj, 2),1)) / signalSize;

stats = extractStatistics(meanIntensity, settings);

end

