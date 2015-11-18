settings = settings();
close all
clear all

% discs = [1:5, 12:17, 21, 22, 24:31];
% AP = {'A', 'P'};
% discs=2;
% AP = {'A'};

filename = 'rawDiscData2/Raw/002.tif';
filename = 'rawDiscData2/Segmented/002P.tif';

% Read video
imageInfo = imfinfo(filename);
numFrames = length(imageInfo);
imSize = [imageInfo(1).Height,imageInfo(1).Width,numFrames];
video = zeros(imSize);
for frame = 1:numFrames
    video(:,:,frame) = imread(filename,frame);
end

% for i = discs
%     for j = 1:length(AP)
%         % Read video
%         filename = ['rawDiscData2/Segmented/'  num2str(i,'%03d') AP{j} '.tif'];
%         imageInfo = imfinfo(filename);
%         numFrames = length(imageInfo);
%         imSize = [imageInfo(1).Height,imageInfo(1).Width,numFrames];
%         video = zeros(imSize);
%         for frame = 1:numFrames
%             video(:,:,frame) = imread(filename,frame);
%         end
%     end 
% end

j = 1;
zscale = 2;
xscale = 50;
yscale = 50;

for i = 1:zscale:(361-zscale)
    rawData(:,:,j) = imresize(mean(video(:,:,i:(i+zscale-1)),3),[round(512/xscale),round(512/yscale)]);
    j = j + 1;
end

a = rawData;
b = convn(a,a(end:-1:1,end:-1:1,end:-1:1));

convSize = size(b);
midTime = round(convSize(3)/2);
convPeak = max(b(:));
%b = squeeze(b(round(convSize(1)/2),round(convSize(2)/2),:))';
b = squeeze(max(max(b,[],2),[],1))';
b(1:midTime) = b(1:midTime) ./ (1:midTime)*(convPeak/midTime);
b(midTime+1:end) = b(midTime+1:end) ./ fliplr(1:(convSize(3)-midTime))*(convPeak/midTime);
plot(b);
b = b(round(midTime/8):round(15*midTime/8));
[pks,locs,w,p] = findpeaks(b);
sort(p)
findpeaks(b);

locations = locs(p>(0.2*max(p(:))));
distances = [locations 0] - [0 locations];
median(distances(2:(end-1)))
% handle = implay(b);
% handle.Visual.ColorMap.UserRangeMin = min(b(:));
% handle.Visual.ColorMap.UserRangeMax = max(b(:));
% handle.Visual.ColorMap.UserRange = 1;
% handle.Visual.ColorMap.MapExpression = 'jet';

dim = (size(rawData));
v = VideoWriter('yesSubtraction.avi');
open(v);

minData = min(rawData,[],3);

for t = 1:dim(3)
    subData(:,:,t) = rawData(:,:,t) - minData;
end

map = 'jet';
n = 3;                   % Decomposition Level
w = 'sym4';              % Near symmetric wavelet
WT = wavedec3(subData,n,w);    % Multilevel 3D wavelet decomposition.

A = cell(1,n);
D = cell(1,n);
for k = 1:n
    A{k} = waverec3(WT,'a',k);   % Approximations (low-pass components)
    D{k} = waverec3(WT,'d',k);   % Details (high-pass components)
end

idxImages_New = [1:dim(3)];
nbIMG = length(idxImages_New);
for ik = 1:nbIMG
    j = idxImages_New(ik);
%     figure('DefaultAxesXTick',[],'DefaultAxesYTick',[],...
%         'DefaultAxesFontSize',8,'Color','w')
    
    temp = [];
    for k = 1:n
        image1 = A{k}(:,:,j);
        image2 = abs(D{k}(:,:,j));
        image2 = (image2 - min(image2(:)))/(max(image2(:)) - min(image2(:)))*max(image1(:));
        temp2 = cat(1,image1,image2);
        temp = cat(2,temp,temp2);
    end
    temp2 = cat(1,rawData(:,:,ik)/max(max(rawData(:,:,ik)))*max(max(subData(:,:,ik))),subData(:,:,ik));
    temp = cat(2,temp,temp2);
    imshow(temp,[])
    colormap(map)
    drawnow
    writeVideo(v,getframe);
end
close(v)










% output = cat(2, video);
% implay(output)
        
