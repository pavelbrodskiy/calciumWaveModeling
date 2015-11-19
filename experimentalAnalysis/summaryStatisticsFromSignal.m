clear all
close all

addpath('downloadedDependencies','customDependencies');
settings = settings();

discs = [1:5, 12:17, 21, 22, 24:31];
AP = {'A', 'P'};

% Stats is a 3D matrix where i represents the wing disc number, j
% represents the compartment, and k represents the bin within the
% compartment.
mmmm1 = 0; nnnn1 = 0;
for mmmm = 0.5:-0.1:0
    for nnnn = [5,1,2,10,25,3,50,70,1000,100]
        nnnn1 = nnnn1 + 1;
        tic
        settings.cutoff=mmmm;
        settings.binDimension=nnnn;
        mmmm1 = mmmm1 + 1;
        
        ii = 1;
        discsDone = 0;
        for i = discs
            for j = 1:length(AP)
                disp(['Processing disc number ' num2str(i,'%03d') ', ' AP{j} ' compartment. ' num2str(100*discsDone/(length(discs)*length(AP)),'%.2f') '% done.'])
                
                % Read video
                filename = ['experimentalData/Segmented/'  num2str(i,'%03d') AP{j} '.tif'];
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
                    temp2 = extractStatistics(videoBins{k}, settings);
                    %settings.outputFunction([num2str(i) ' ' num2str(j) ' ' num2str(k)]);
                    if temp2.flag
                        statTemp = [statTemp temp2];
                    end
                end
                stats{ii,j} = statTemp;
                discsDone = discsDone + 1;
            end
            ii = ii + 1;
        end
        
        analysis = analyzeStats(stats, settings);
        close all
        
        if ~isempty(analysis)
            p1111(mmmm1,nnnn1) = analysis.pFrequency
            t1111(mmmm1,nnnn1) = toc;
            s1111{mmmm1,nnnn1} = stats;
            n1111(mmmm1,nnnn1) = analysis.n;
            save('dataDump',analysis.pFrequency,toc,stats,analysis.n,nnnn,mmmm,'-append')
        end
    end
end