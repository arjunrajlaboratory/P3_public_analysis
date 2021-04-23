function [centroidsAndTheirSpotCounts] = extractSpotsAndCentroids(dataTopDir, wellPath)
% wellPath is a child (or grandchild) of dataTopDir/

cd([dataTopDir, wellPath])

spotsAndCentroids = dentist.load();

centroidLocations = spotsAndCentroids.getCentroids;
alexanum = spotsAndCentroids.getNumSpotsForCentroids('alexa');
cynum = spotsAndCentroids.getNumSpotsForCentroids('cy');
tmrnum = spotsAndCentroids.getNumSpotsForCentroids('tmr');

centroidsAndTheirSpotCounts = table(centroidLocations.xPositions, ...
    centroidLocations.yPositions, ...
    alexanum, ...
    tmrnum, ...
    cynum, ...
    [1:size(tmrnum,1)]');
centroidsAndTheirSpotCounts.Properties.VariableNames = {'xPos', 'yPos', 'alexa', 'tmr', 'cy', 'cellID'};

end