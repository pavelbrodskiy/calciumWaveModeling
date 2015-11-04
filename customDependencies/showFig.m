function showFig( concentration, t )

concentration(isnan(concentration)) = 0;

imshow(imresize(concentration,[1000,1000]),[]);
colormap('jet');

%imwrite(uint16(concentration*1000), [t ' data.tif']);
concentration = imresize(concentration,[2000,2000]);
concentration = concentration - min(concentration(:));
concentration = concentration ./ max(concentration(:));
concentration = concentration * 256;
imwrite(grs2rgb(concentration,jet(256)) , [t ' image.png']);

end

