plot(t,meanIntensity{4,2})

title(['Wing Disc ' num2str(4, '%03d')])
    ylabel('Intensity');
    xlabel('Time (min)');
    
    
    print(['Wing Disc ' num2str(4, '%03d')],'-painters', '-dpng', '-r1200')