function showFig( concentration, t, p )

concentration(isnan(concentration)) = 0;

concentration = imresize(concentration,[p.kymDim, p.kymDim]);

imwrite(uint16(concentration*10000), [t 'rawData.png']);
concentration = concentration - min(concentration(:));
concentration = concentration ./ max(concentration(:));
concentration = concentration * 256;
image = grs2rgb(flipud(concentration),jet(256));
imwrite(image, [t 'image.png']);

if any(p.outputModes==5)
    imshow(image,[]);
end

end

