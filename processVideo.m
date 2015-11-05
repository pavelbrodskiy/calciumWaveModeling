function [ analysisObj ] = processVideo( varargin )

filename = 'WingDisc2 P compartment.tif';

if nargin > 0
    video = varargin{1};
else
    imageInfo = imfinfo(filename);
    numFrames = length(imageInfo);
    imSize = [imageInfo(1).Height,imageInfo(1).Width,numFrames];
    video = zeros(imSize);
    for frame = 1:numFrames
        video(:,:,frame) = imread(filename,frame);
    end
end


for i = 1:numFrames
    subtractedVideo(:,:,i) = video(:,:,i) - min(video,[],3);
    signal = video(:,:,i);
    %cleanedSignal = subtractedVideo(:,:,i);
    analysisObj.meanIntensity(i) = mean(signal(:));
    %cleanedIntensity(i) = mean(cleanedSignal(:));
    analysisObj.stdevIntensity(i) = std(signal(:));
end



end

