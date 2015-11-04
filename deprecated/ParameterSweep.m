% Add '$P' to front of parameter if only P compartment

runNumber = 200;
percentVariation = 1;
timeAtStart = now();

sweeper( 1.5, 'noiseSigma', percentVariation, runNumber, timeAtStart);

