statData = stats;

AFrequency = [statData(:,1).meanFrequency]./60;
PFrequency = [statData(:,2).meanFrequency]./60;
AWidth = [statData(:,1).meanWidth]*60;
PWidth = [statData(:,2).meanWidth]*60;
AAmplitude = [statData(:,1).meanAmplitude];
PAmplitude = [statData(:,2).meanAmplitude];

[~,pFrequency,~,~] = ttest(AFrequency, PFrequency);
[~,pWidth,~,~] = ttest(AWidth, PWidth);
[~,pAmplitude,~,~] = ttest(AAmplitude, PAmplitude);

labels = [repmat({'A Compartment'},1,21), repmat({'P Compartment'},1,21)];

%boxplot([AFrequency, PFrequency], labels);
close all
notBoxPlot([AFrequency; PFrequency]',[1,1.5])
title('Frequency A vs P')
ylabel('Frequency (Hz)');
xlim([0.8, 1.7])
print('Frequency A vs P', '-painters', '-dpng', '-r1200')

%boxplot([AWidth, PWidth], labels);
close all
notBoxPlot([AWidth; PWidth]',[1,1.5])
title('Width A vs P')
ylabel('Width (s)');
xlim([0.8, 1.7])
print('Width A vs P', '-painters', '-dpng', '-r1200')

%boxplot([AAmplitude, PAmplitude], labels);
close all
notBoxPlot([AAmplitude; PAmplitude]',[1,1.5])
title('Amplitude A vs P')
ylabel('Amplitude (AU)');
xlim([0.8, 1.7])
print('Amplitude A vs P', '-painters', '-dpng', '-r1200')
