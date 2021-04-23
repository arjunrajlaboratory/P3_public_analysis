function [] = extractStackSpotCountsAndMasks(dataTopDir, wellPath, outDir)
% wellPath is a child (or grandchild) of dataTopDir/

cd([dataTopDir, wellPath])

dataExtractor = improc2.launchDataExtractor;
 
dataExtractor.extractFromProcessorData('cy.RNACounts', @getNumSpots, 'cy')
dataExtractor.extractFromProcessorData('alexa.RNACounts', @getNumSpots, 'alexa')
dataExtractor.extractFromProcessorData('tmr.RNACounts', @getNumSpots, 'tmr')

dataExtractor.extractAllToCSVFile([outDir, '/counts.csv'])

tools = improc2.launchImageObjectTools;

objectNumber = [];
arrayNumber = [];
isGood = [];
nPixels = [];
x = table(objectNumber, arrayNumber, isGood, nPixels);

while tools.iterator.continueIteration
    
    objectNumber = tools.navigator.currentObjNum();
    arrayNumber = tools.navigator.currentArrayNum();
    mask = tools.objectHandle.getCroppedMask;
    nPixels = sum(sum(mask));
    isGood = tools.annotations.getValue('isGood');
    
    x = [x; table(objectNumber, arrayNumber, isGood, nPixels)];
    
    tools.iterator.goToNextObject()
end

writetable(x, [outDir, '/cellMaskSize.csv']);
end