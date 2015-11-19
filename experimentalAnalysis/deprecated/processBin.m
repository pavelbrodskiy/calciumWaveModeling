function stats = processBin(videoBin, settings)
% This function takes a video bin and converts it into a 1D signal, then
% extracts summary statistics.

if nargin < 2
    settings = settings();
end

resizedVideo = videoResize(videoBin, settings.rescaleBy);

stats = extractStatistics(meanIntensity, settings);

end

