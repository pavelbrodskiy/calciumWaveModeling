function videoBins = videoCutter(video, settings)

[x, y, ~] = size(video);

maxProj = max(video,[],3);
regionShape = squeeze(maxProj > 0);

binNumber = 1;
for i = 1:settings.binDimension:(x-settings.binDimension)
    for j = 1:settings.binDimension:(y-settings.binDimension)
        if min(min(regionShape(i:(i+settings.binDimension-1),(j:(j+settings.binDimension-1)))))
            videoBins{binNumber} = video(i:(i+settings.binDimension-1),(j:(j+settings.binDimension-1)),:);
            binNumber = binNumber + 1;
        end
    end
end

if binNumber == 1
    videoBins{1} = video;
end

end