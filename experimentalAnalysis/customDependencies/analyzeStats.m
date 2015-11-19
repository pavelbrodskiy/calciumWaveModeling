function analysis = analyzeStats(stats, settings)


%% Extract statistics
[x, ~] = size(stats);

analysis.n = 0;

statData = processCellMatrix(stats);
% a = isnan(statData);
if isempty(statData)
    analysis.pFrequency = 1;
    analysis.pWidth = 1;
    analysis.pAmplitude = 1;
else
    b = cellfun(@isempty,statData)
    c = b(:,1)|b(:,2);
    statData(c,:) = [];
end

if isempty(statData)
    analysis.pFrequency = 1;
    analysis.pWidth = 1;
    analysis.pAmplitude = 1;
else
    
    statMat = cell2mat(statData);
    
    AFrequency = [statMat(:,1).meanFrequency]./60;
    PFrequency = [statMat(:,2).meanFrequency]./60;
    AWidth = [statMat(:,1).meanWidth]*60;
    PWidth = [statMat(:,2).meanWidth]*60;
    AAmplitude = [statMat(:,1).meanAmplitude];
    PAmplitude = [statMat(:,2).meanAmplitude];
    
    [~,analysis.pFrequency,~,~] = ttest(AFrequency, PFrequency);
    [~,analysis.pWidth,~,~] = ttest(AWidth, PWidth);
    [~,analysis.pAmplitude,~,~] = ttest(AAmplitude, PAmplitude);
    [analysis.n, ~] = size(statMat);
    
    %% Output if needed
    labels = [repmat({'A Compartment'},1,analysis.n), repmat({'P Compartment'},1,analysis.n)];
    
    if(settings.outputBoxPlot)
        close all
        boxplot([AFrequency, PFrequency], labels);
        title(['Frequency A vs P    p = ' num2str(analysis.pFrequency)])
        ylabel('Frequency (Hz)');
        print([settings.outputDirectory 'Frequency A vs P box'], '-painters', '-dpng', '-r1200')
        
        close all
        boxplot([AWidth, PWidth], labels);
        title(['Width A vs P    p = ' num2str(analysis.pWidth)])
        ylabel('Width (s)');
        print([settings.outputDirectory 'Width A vs P box'], '-painters', '-dpng', '-r1200')
        
        close all
        boxplot([AAmplitude, PAmplitude], labels);
        title(['Amplitude A vs P    p = ' num2str(analysis.pAmplitude)])
        ylabel('Amplitude (AU)');
        print([settings.outputDirectory 'Amplitude A vs P box'], '-painters', '-dpng', '-r1200')
        
    end
    
    if(settings.outputNotBoxPlot)
        close all
        notBoxPlot([AFrequency; PFrequency]',[1,1.5])
        title(['Frequency A vs P    p = ' num2str(analysis.pFrequency)])
        ylabel('Frequency (Hz)');
        xlim([0.8, 1.7])
        print([settings.outputDirectory 'Frequency A vs P'], '-painters', '-dpng', '-r1200')
        
        close all
        notBoxPlot([AWidth; PWidth]',[1,1.5])
        title(['Width A vs P    p = ' num2str(analysis.pWidth)])
        ylabel('Width (s)');
        xlim([0.8, 1.7])
        print([settings.outputDirectory 'Width A vs P'], '-painters', '-dpng', '-r1200')
        
        close all
        notBoxPlot([AAmplitude; PAmplitude]',[1,1.5])
        title(['Amplitude A vs P    p = ' num2str(analysis.pAmplitude)])
        ylabel('Amplitude (AU)');
        xlim([0.8, 1.7])
        print([settings.outputDirectory 'Amplitude A vs P'], '-painters', '-dpng', '-r1200')
    end
    
    if(settings.outputRatioPlot)
        boxplot([AFrequency, PFrequency], labels);
        close all
        notBoxPlot([(AFrequency./PFrequency); (AWidth./PWidth); (AAmplitude./PAmplitude)]',[1,1.5,2])
        title('Frequency, Width, Amplitude Ratios')
        ylabel('Ratio');
        xlim([0.8, 2.2, ])
        print([settings.outputDirectory 'Frequency, Width, Amplitude Ratios'], '-painters', '-dpng', '-r1200')
    end
    
    if(settings.outputScatter)
        close all
        
        subplot(2,2,1)
        scatter(AFrequency,PFrequency)
        title('Frequency')
        xlabel('Anterior')
        ylabel('Posterior')
        
        subplot(2,2,2)
        scatter(AWidth,PWidth)
        title('Width')
        xlabel('Anterior')
        ylabel('Posterior')
        
        subplot(2,2,3)
        scatter(AAmplitude,PAmplitude)
        title('Amplitude')
        xlabel('Anterior')
        ylabel('Posterior')
        
        print([settings.outputDirectory 'Frequency, Width, Amplitude Scatter'], '-painters', '-dpng', '-r1200')
    end
end
end

