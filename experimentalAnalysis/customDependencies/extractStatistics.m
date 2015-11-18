function stats = extractStatistics(signalOverTime, settings)
% This function takes a 1D signal and extracts summary statistics.

if nargin < 2
    settings = settings();
end

[acf,lags,bounds] = autocorr(signalOverTime,360);
[~,peakTimes,widths,amplitudes] = findpeaks(acf);

amptemp = sort(amplitudes,'descend');
% amptemp
% amplitudes == amptemp(2)
% peakTimes(amplitudes == amptemp(2))
lagInSeconds = peakTimes(amplitudes == amptemp(2));
stats.meanFrequency = 1/lagInSeconds;

peakTimes = [1, 2];

% levels = statelevels(signalOverTime);
% [Rise,LoTime,HiTime,LoLev,HiLev] = risetime(signalOverTime,settings.t);
% [~,peakTimes,widths,amplitudes] = findpeaks(signalOverTime);
% tempAmplitudes = sort(amplitudes, 'descend');
% cutoff = tempAmplitudes(2) * settings.amplitudeThreshold;
% 
% peakTimes(amplitudes < cutoff) = [];
% widths(amplitudes < cutoff) = [];
% amplitudes(amplitudes < cutoff) = [];

if length(peakTimes) > 1
    timeBetweenPeaks = peakTimes(2:end) - peakTimes(1:end-1);
%     stats.meanFrequency = mean(1./timeBetweenPeaks);
    stats.stdevFrequency = std(1./timeBetweenPeaks);
    stats.meanAmplitude = mean(amplitudes);
    stats.stdevAmplitude = std(amplitudes);
    stats.meanWidth = mean(widths);
    stats.stdevWidth = std(widths);
    stats.stdevFrequency = 1;
    stats.meanAmplitude = 1;
    stats.stdevAmplitude = 1;
    stats.meanWidth = 1;
    stats.stdevWidth = 1;
    stats.flag = true;
else
    stats.flag = false;
end


end