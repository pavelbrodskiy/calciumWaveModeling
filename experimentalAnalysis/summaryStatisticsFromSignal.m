clear all
close all

settings = settings();

discs = [1:5, 12:17, 21, 22, 24:31];
AP = {'A', 'P'};
addpath('downloadedDependencies','customDependencies');

% Stats is a 3D matrix where i represents the wing disc number, j
% represents the compartment, and k represents the bin within the
% compartment.

ii = 1;
for i = discs
    for j = 1:length(AP)
        [ii, j]
        % Read video
        filename = ['rawDiscData2/Segmented/'  num2str(i,'%03d') AP{j} '.tif'];
        imageInfo = imfinfo(filename);
        numFrames = length(imageInfo);
        imSize = [imageInfo(1).Height,imageInfo(1).Width,numFrames];
        video = zeros(imSize);
        for frame = 1:numFrames
            video(:,:,frame) = imread(filename,frame);
        end
        
        % Bin video
        videoBins = videoCutter(video, settings);
        
        % Analyze video bins
        statTemp = [];
        for k = 1:length(videoBins)
            if temp2.flag
                statTemp = [statTemp extractStatistics(videoBins{k}, settings)];
            end
        end
        stats{ii,j} = statTemp;
    end
    ii = ii + 1;
end

analysis = analyzeStats(stats, settings);
close all

