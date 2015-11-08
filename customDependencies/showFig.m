function showFig( concentration, t, p )

concentration(isnan(concentration)) = 0;

concentration = imresize(concentration,p.kymDim);
concentration = concentration - min(concentration(:));
concentration = concentration ./ max(concentration(:));
concentration = concentration * 256;
imwrite(grs2rgb(flipud(concentration), jet(256)) , [t 'image.png']);

end

