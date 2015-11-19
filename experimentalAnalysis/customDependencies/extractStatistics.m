function stats = extractStatistics(compressedVideo, settings)
% This function takes a 3D signal and extracts summary statistics.

if nargin < 2
    settings = settings();
end

maxProj = max(compressedVideo,[],3);
minProj = mean(mean(min(compressedVideo,[],3)));
signalSize = sum(sum(maxProj > 0));
meanIntensity = squeeze(sum(sum(compressedVideo - minProj, 2),1)) / signalSize;

correlation = autocorr(meanIntensity,length(meanIntensity)-1);

[~,locs,~,p] = findpeaks(correlation);
[~,locs2,~,p2] = findpeaks(meanIntensity);

peakHeights = sort(p,'descend');
peakHeights2 = sort(p2,'descend');

timeBetweenPeaks = [locs(p == peakHeights(1)) locs(p == peakHeights(1))];
xs = 1:length(correlation);
ys = xs*0;
firstPeak = locs2(p2==peakHeights2(1));
secondPeak = firstPeak + timeBetweenPeaks(1);
thirdPeak = firstPeak - timeBetweenPeaks(1);
thirdPeak = max([1, thirdPeak]);
ys([firstPeak, secondPeak, thirdPeak]) = max(meanIntensity(:));
ys = ys(1:length(correlation));
plot(xs,meanIntensity,xs,ys)
drawnow
disp(num2str(timeBetweenPeaks));

if length(timeBetweenPeaks) > 1
    stats.meanFrequency = mean(1./timeBetweenPeaks);
    stats.stdevFrequency = std(1./timeBetweenPeaks);
    stats.meanAmplitude = [1 1];
    stats.stdevAmplitude = [1 1];
    stats.meanWidth = [1 1];
    stats.stdevWidth = [1 1];

    stats.flag = true;
else
    stats.flag = false;
end

end