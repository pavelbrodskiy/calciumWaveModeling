function smallVideo = videoResize(rawVideo, scale)

% Input a 3D matrix and scaling factors. This function outputs the video
% with smaller dimensions scaled by the scaling factors.

j = 1;
for i = 1:scale(3):(361-scale(3))
    smallVideo(:,:,j) = imresize(mean(rawVideo(:,:,i:(i+scale(3)-1)),3),[round(512/scale(1)),round(512/scale(2))]);
    j = j + 1;
end

end