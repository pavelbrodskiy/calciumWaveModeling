clear all
close all

discs = [1:5, 12:17, 21, 22, 24:31];


m = 1
for k = discs
    j = 1;
    for tmpStr = {'A', 'P'}
        filename = ['rawDiscData2/Segmented/'  num2str(k,'%03d') tmpStr{1} '.tif'];
        imageInfo = imfinfo(filename);
        numFrames = length(imageInfo);
        imSize = [imageInfo(1).Height,imageInfo(1).Width,numFrames];
        video = zeros(imSize);
        for frame = 1:numFrames
            video(:,:,frame) = imread(filename,frame);
        end
        
        maxProj = max(video,[],3);
        signalSize = sum(sum(maxProj > 0));
        meanIntensity{m,j} = squeeze(sum(sum(video,2),1)) / signalSize;
        
        j = j + 1;
    end
    m = m + 1
end

t = 0:0.1666:60;

% for i = 1:6
% subplot(3,2,i);
% plot(t,meanIntensity{i,1},t,meanIntensity{i,2});
% title(['Wing Disc ' num2str(i)])
% ylabel('Average Intensity (AU)');
% xlabel('Time (min)');
% end

% for i = 1:(m-1)
%     subplot(5,5,i);
%     plot(t,(meanIntensity{i,1}-min(meanIntensity{i,1}))/(max(meanIntensity{i,1})-min(meanIntensity{i,1})),t,(meanIntensity{i,2}-min(meanIntensity{i,2}))/(max(meanIntensity{i,2})-min(meanIntensity{i,2})));
%     title(['Wing Disc ' num2str(discs(i), '%03d')])
%     ylabel('Normalized Intensity');
%     xlabel('Time (min)');
% end
% 
% subplot(5,5,m);
% plot(1,1,1,1);
% legend('A Compartment','P Compartment','Location','best')
% 
% 
% %print('Wing Disc Compartment Analysis','-depsc')
% print('Wing Disc Compartment Analysis 2','-painters', '-dpdf', '-r1200')
if(false)
for i = 1:(m-1)
   
    plot(t,meanIntensity{i,1},t,meanIntensity{i,2});
    title(['Wing Disc ' num2str(discs(i), '%03d')])
    ylabel('Intensity');
    xlabel('Time (min)');
    
    legend('A Compartment','P Compartment','Location','best')
    
    print(['Wing Disc ' num2str(discs(i), '%03d')],'-painters', '-dpdf', '-r1200')
end

for i = 1:(m-1)
   
    plot(t,(meanIntensity{i,1}-min(meanIntensity{i,1}))/(max(meanIntensity{i,1})-min(meanIntensity{i,1})),t,(meanIntensity{i,2}-min(meanIntensity{i,2}))/(max(meanIntensity{i,2})-min(meanIntensity{i,2})));
    title(['Wing Disc ' num2str(discs(i), '%03d')])
    ylabel('Normalized Intensity');
    xlabel('Time (min)');
    
    legend('A Compartment','P Compartment','Location','best')
    
    print(['Normalized Wing Disc ' num2str(discs(i), '%03d')],'-painters', '-dpdf', '-r1200')
end

end

%print('Wing Disc Compartment Analysis','-depsc')


