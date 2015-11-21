function settings = settings2()

settings.amplitudeThreshold = 0.3;
settings.t = 0:0.1666:60;

settings.outputScatter = false;
settings.outputRatioPlot = false;
settings.outputNotBoxPlot = false;
settings.outputBoxPlot = false;

% Directory to save output graphs to.
settings.outputDirectory = ['analysisOutput' filesep];

% Pixel dimensions of bins. Use dim > x for whole video.
settings.binDimension = 25;

% Ratio to divide x, y, and t dimensions by.
settings.rescaleBy = [1,1,1]; 

settings.outputFunction = @(outputName) print([settings.outputDirectory outputName], '-painters', '-dpng', '-r10');

settings.cutoff = 0.3;

end