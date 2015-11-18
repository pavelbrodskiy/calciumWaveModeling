function settings = settings()

settings.amplitudeThreshold = 0.3;
settings.t = 0:0.1666:60;

settings.outputScatter = false;
settings.outputRatioPlot = true;
settings.outputNotBoxPlot = true;
settings.outputBoxPlot = false;

% Directory to save output graphs to.
settings.outputDirectory = ['analysisOutput' filesep];

% Pixel dimensions of bins. Use dim > x for whole video.
settings.binDimension = 1055;

% Ratio to divide x, y, and t dimensions by.
settings.rescaleBy = [50,50,10]; 

end