%% Clear the interface
clc;
clear all;
close all;

%% Parameters
discs = fliplr([1:5, 12:17, 21, 22, 24:31]);
mmm = 1;
K=8;
    
for kkk = discs
    %% Read in tiff stack
    filename = ['rawDiscData2/Raw/'  num2str(kkk,'%03d') '.tif'];
    imageInfo = imfinfo(filename);
    numFrames = length(imageInfo);
    imSize = [imageInfo(1).Height,imageInfo(1).Width,numFrames];
    X = zeros(imSize);
    for frame = 1:numFrames
        X(:,:,frame) = imread(filename,frame);
    end
    
    %% Subtract mean image for zero-mean deviation
    [h,w,nFrames]=size(X);
    Xmean=mean(X,3);
    for n=nFrames:-1:1
        Ws(:,:,n)=X(:,:,n)-Xmean;
    end
    
    %% Resize via low-pass filtering for SVD analysis
    A=zeros(nFrames,h*w);
    for n=1:nFrames
        x=Ws(:,:,n);
        A(n,:)=x(:)';
    end
    
    [u,s,v]=svds(A,K); % s is the singular value matrix, v is the eigenmode?
    
    V=zeros(h,w,nFrames);
    U=zeros(nFrames,K);
    for k=1:K
        pm=sign(mean(v(:,k)));
        V(:,:,k)=reshape(pm*v(:,k),h,w);
        U(:,k)=pm*u(:,k);
    end
    
    close all
    figure;
    for k=1:K
        subplot(2,4,k);
        imshow(V(:,:,k),[min(V(:)) max(V(:))]);
        title(sprintf('k=%d',k));
    end
    
    print(['Eigenmodes ' num2str(kkk, '%03d')],'-painters', '-dpdf', '-r1200')
    
end

