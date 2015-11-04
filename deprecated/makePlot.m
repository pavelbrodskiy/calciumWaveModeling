%function makePlot( dataMatrix, v7, P_IP3, v8 )

dataMatrix = width;

meanVal = squeeze(nanmean(dataMatrix,1));
stdVal = squeeze(nanstd(dataMatrix,1));

subplot(3,1,1)
errorbar(v7, meanVal(1,1:length(v7)), stdVal(1,1:length(v7)))
subplot(3,1,2)
errorbar(P_IP3, meanVal(2,1:length(P_IP3)), stdVal(2,1:length(P_IP3)))
subplot(3,1,3)
errorbar(v8, meanVal(3,1:length(v8)), stdVal(3,1:length(v8)))


