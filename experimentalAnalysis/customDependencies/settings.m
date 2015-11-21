function settings = settings()

settings.amplitudeThreshold = 0.3;
settings.t = 0:0.1666:60;

settings.outputScatter = true;
settings.outputRatioPlot = true;
settings.outputNotBoxPlot = true;
settings.outputBoxPlot = true;

% Directory to save output graphs to.
settings.outputDirectory = ['analysisOutput' filesep];

% Pixel dimensions of bins. Use dim > x for whole video.
settings.binDimension = 25;

% Ratio to divide x, y, and t dimensions by.
settings.rescaleBy = [1,1,1]; 

settings.outputFunction = @(outputName) print([settings.outputDirectory outputName], '-painters', '-dpng', '-r10');

settings.cutoff = 0.3;

end