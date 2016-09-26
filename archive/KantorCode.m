%% Clear the interface
clc;
clear all;
close all;

%% Read and convert .avi file
mov = VideoReader('20 percent FEX, FEX Gradient test_s1_t1.TIF - Stage3 "D3"-1.avi');

nFrames=mov.NumberOfFrames;
h=mov.Height;
w=mov.Width;
X(1:nFrames)=struct('cdata',zeros(h,w,3,'uint8'),'colormap',[]);
for n=1:nFrames
    X(n).cdata=read(mov,n);
end

%% Crop by inspection
hcrop=1:512;
wcrop=1:512;

%% Crop and extract luminosity
% figure(2);clf;
X=zeros(length(hcrop),length(wcrop),nFrames);

for n=1:nFrames
    x=double(read(mov,n))/225;
    x=x(hcrop,wcrop,:);
    X(:,:,n)=0.21*x(:,:,1)+0.71*x(:,:,2)+0.08*x(:,:,3);
    %Show the images to verify
%     imshow(X(:,:,n));
%     drawnow;
end

% Generate background image as the minimum luminosity
Xdark=min(X,[],3);
imshow(Xdark);

% Show waves with subtracted background
Y=zeros(size(X));
figure(2);clf;
for n=1:nFrames
    Y(:,:,n)=X(:,:,n)-Xdark;
%     imshow(Y(:,:,n));
%     title(sprintf('n=%d',n));
%     drawnow;
end

%% Create isosurface rendering
Z=Y.*(Y>=0.25); %???
sf=1/4;
Zs=zeros(512/4,512/4,nFrames);
for n=1:nFrames
    Zs(:,:,n)=imresize(Z(:,:,n),sf);
end
XX=[1:1:512/4]*1.37;
YY=XX;
ZZ=[1:1:nFrames]/6;

figure(1);clf;
isosurface(XX,YY,ZZ,Zs,0.25);  %???
xlabel('Vertical/µm');
xlim([0 512*0.3425]);
ylim([0 512*0.3425]);
zlim([0 nFrames/6]);
ylabel('Horizontal/µm');
zlabel('Time/min');
axis vis3d
grid

% 3D Movie
% OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
% CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'WellMadeVid',OptionZ)

%% Subtract mean image for zero-mean deviation
[h,w,nFrames]=size(Y);
W=zeros(h,w,nFrames);
Ymean=mean(Y,3);
for n=1:nFrames
    W(:,:,n)=Y(:,:,n)-Ymean;
end

Wmin=min(W(:));
Wmax=max(W(:));

figure(2);
for n=1:nFrames
    imshow(W(:,:,n)+abs(Wmin)); %???? Why add the minimum value to every element
%     title(sprintf('n=%d',n))
%     drawnow;
end

%% Image Covariances
% C=zeros(nFrames,nFrames);
% for i=1:nFrames
%     for j=i:nFrames
%         c=W(:,:,i).*W(:,:,j);
%         C(i,j)=sum(c(:));
%         C(j,i)=C(i,j);
%     end
% end
% 
% contour(C);
% axis('equal');
% xlabel('Image Number');
% ylabel('Image Number');

%% Resize via low-pass filtering for SVD analysis
sf=1/2;
Ws=zeros(h/2,w/2,nFrames);
for n=1:nFrames
    Ws(:,:,n)=imresize(W(:,:,n),sf);
end

figure(3);
for n=1:nFrames
    imshow(Ws(:,:,n)+abs(Wmin));
    title(sprintf('n=%d',n));
    drawnow;
end


K=8;

[h,w,nFrames]=size(Ws);

A=zeros(nFrames,h*w);
for n=1:nFrames
    x=Ws(:,:,n);
    A(n,:)=x(:)';  %??? What does ' mean?
end

[u,s,v]=svds(A,K); % s is the singular value matrix, v is the eigenmode?

V=zeros(h,w,nFrames);
U=zeros(nFrames,K);
for k=1:K
    pm=sign(mean(v(:,k)));
    V(:,:,k)=reshape(pm*v(:,k),h,w);
    U(:,k)=pm*u(:,k);
end

figure(4);
subplot(2,4,1);
for k=1:K
    subplot(2,4,k);
    imshow(0.5+20*V(:,:,k));
    title(sprintf('k=%d',k));
end


% %Open output video
% figure(4);clf;
% writeObj=VideoWriter('CalciumwaveAnalysis.avi')
% open(writeObj);

%Construct trajectory
g=U*s;

% %Plot trajectory
% subplot(2,2,[3 4]);
% plot(g)
% hold on
% ax=axis;
% h1=plot([0 0],[ax(3) ax(4)],'r');
% hold off
% xlabel('Frame Number');
% title(sprintf('%d mode reconstruction',K));
% lgd={};
% for k=1:K
%     lgd{k}=sprintf('k=%d',k);
% end
% legend(lgd,'location','EastOutside');

% subplot(2,4,1);
% for n=1:nFrames
%     subplot(2,2,1);
%     imshow(X(:,:,n));
%     subplot(2,2,2);
%     Wa=zeros(h,w);
%     for k=1:K
%         Wa=Wa+g(n,k)*V(:,:,k);
%     end
%     Ya=Ymean+imresize(Wa,1/sf);
%     Xa=Xdark+Ya;
%     imshow(Xa);
%     title(sprintf('Reconstruction Frame=%d',n))
%     
%     subplot(2,2,[3 4]);
%     set(h1,'XData',[n n]);
%     
%     frame=getframe(gcf);
%     writeVideo(writeObj,frame);
%     drawnow
% end
% 
% close(writeObj);

%Phase plot for eigenmodes
figure(5);
plot3(g(:,1),g(:,2),g(:,3),'Linewidth',2);
grid
xlabel('Mode 1');
ylabel('Mode 2');
zlabel('Mode 3');

%% Autocorrelation
% figure(6);
% for k=1:min([K,4])
%     subplot(min([K,4]),1,k);
%     autocorr(g(:,k),300);
%     title(sprintf('Autocorrelation mode k=%d',k));
% end
% 
% %% Crosscorrelation
% figure(7);
% subplot(3,1,1);
% crosscorr(g(:,1),g(:,2),300)
% title('Cross Correlation 1 vs 2');
% 
% subplot(3,1,2);
% crosscorr(g(:,1),g(:,3),300)
% title('Cross Correlation 1 vs 3');
% 
% subplot(3,1,3);
% crosscorr(g(:,2),g(:,3),300)
% title('Cross Correlation 2 vs 3')