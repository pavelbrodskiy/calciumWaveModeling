subplot(3,1,1)
plot(signalOverTime)
subplot(3,1,2)
plot(gradient(signalOverTime))
subplot(3,1,3)
plot(del2(signalOverTime))

% [xcf,lags,bounds] = crosscorr(signalOverTime
[acf,lags,bounds] = autocorr(signalOverTime,360);

plot(acf)

% 
% video1=video(:,:,1:100);
% video2=video(:,:,101:200);
% video3=video(:,:,201:300);
% video4=video(:,:,301:end);
% 
% video1(video1==0)=[];
% video2(video2==0)=[];
% video3(video3==0)=[];
% video4(video4==0)=[];
% 
% hist(video1,100);
% axis([0, 7e3,0, 8e5]);
% hist(video2,100);
% axis([0, 7e3,0, 8e5]);
% hist(video3,100);
% axis([0, 7e3,0, 8e5]);
% hist(video4,100);
% axis([0, 7e3,0, 8e5]);
% 0