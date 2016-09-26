for k = randperm(length(discs))
    for comp = randperm(2)
        close all
        plot(t,meanIntensity{k,comp});
        peakTimesManual = ginput();
        peakTimesManual = sort(peakTimesManual(:,1));
        timeBetweenPeaks = peakTimesManual(2:end) - peakTimesManual(1:end-1);
        [~,peakTimes,widths,amplitudes] = findpeaks(meanIntensity{k,comp}, 6, 'MinPeakDistance', 0.9 * min(timeBetweenPeaks)); %'NPeaks', length(peakTimesManual), 
        
        tempAmp = sort(amplitudes, 'descend');
        
        peakTimes(amplitudes < tempAmp(min(length(amplitudes),length(peakTimesManual)))) = [];
        widths(amplitudes < tempAmp(min(length(amplitudes),length(peakTimesManual)))) = [];
        amplitudes(amplitudes < tempAmp(min(length(amplitudes),length(peakTimesManual)))) = [];
        
        if length(peakTimes) > 1
            
            timeBetweenPeaks = peakTimes(2:end) - peakTimes(1:end-1);
            stats(k,comp).meanFrequency = mean(1./timeBetweenPeaks);
            stats(k,comp).stdevFrequency = std(1./timeBetweenPeaks);
            stats(k,comp).meanAmplitude = mean(amplitudes);
            stats(k,comp).stdevAmplitude = std(amplitudes);
            stats(k,comp).meanWidth = mean(widths);
            stats(k,comp).stdevWidth = std(widths);
        
        else
            stats(k,comp) = NaN;
        end
        
        
    end
end