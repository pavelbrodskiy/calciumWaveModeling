%% Clear the interface
clc;
clear all;
close all;

discs = [1:5, 12:17, 21, 22, 24:31];
mmm = 1;

for kkk = discs
    filename = ['rawDiscData2/Raw/'  num2str(kkk,'%03d') '.tif'];
    imageInfo = imfinfo(filename);
    numFrames = length(imageInfo);
    imSize = [imageInfo(1).Height,imageInfo(1).Width,numFrames];
    X = zeros(imSize);
    for frame = 1:numFrames
        X(:,:,frame) = imread(filename,frame);
    end
    
    % Generate background image as the minimum luminosity
    Xdark=min(X,[],3);
    
    Y=X-repmat(Xdark,[1,1,361]);
    
    %% Subtract mean image for zero-mean deviation
    [h,w,nFrames]=size(Y);
    W=zeros(h,w,nFrames);
    Ymean=mean(Y,3);
    for n=1:nFrames
        Ws(:,:,n)=Y(:,:,n)-Ymean;
    end
    
    %% Resize via low-pass filtering for SVD analysis
    
    K=8;
    
    [h,w,nFrames]=size(Ws);
    
    A=zeros(nFrames,h*w);
    for n=1:nFrames
        x=Ws(:,:,n);
        A(n,:)=x(:)';
    end
    
    [u,s,v]=svds(A,K); % s is the singular value matrix, v is the eigenmode?
    
    weightArray = s(:);
    weightArray (weightArray == 0) = [];
    weights{mmm} = weightArray;
    mmm = mmm + 1;
    
    V=zeros(h,w,nFrames);
    U=zeros(nFrames,K);
    for k=1:K
        pm=sign(mean(v(:,k)));
        V(:,:,k)=reshape(pm*v(:,k),h,w);
        U(:,k)=pm*u(:,k);
    end
    
    close all
    figure;
    subplot(2,4,1);
    for k=1:K
        subplot(2,4,k);
        imshow(0.5+20*V(:,:,k));
        title(sprintf('k=%d',k));
    end
    
    print(['Eigenmodes ' num2str(discs(kkk), '%03d')],'-painters', '-dpdf', '-r1200')
    
end

