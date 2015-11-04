function showFig( concentration, t )

concentration(isnan(concentration)) = 0;

concentration = imresize(concentration,[1000,1000]);
concentration = concentration - min(concentration(:));
concentration = concentration ./ max(concentration(:));
concentration = concentration * 256;
imwrite(grs2rgb(flipud(concentration), jet(256)) , [t ' image.png']);

end

