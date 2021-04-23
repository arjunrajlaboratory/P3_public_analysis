function [] = IAM942_34_subset_spotsAndCentroids(dataTopDir, outDir)


experimentDir = 'transdiff/20190731_IAM942-34-subset/';

outExperimentDir = [outDir experimentDir];

if (~exist(outExperimentDir, 'dir'))
    mkdir(outExperimentDir)
end

readme = ReadYaml([dataTopDir, experimentDir, 'readMe.yaml']);

copyfile([dataTopDir, experimentDir, 'readMe.yaml'], [outDir, experimentDir]);

wellMap.wellName = strsplit(readme.wellMap.wellName, ',');
wellMap.condition = strsplit(readme.wellMap.condition, ',');


for i = 1:4
    
    wName = wellMap.wellName{i};
    wellPath = [experimentDir, wName, '_scan'];
    
    results = extractSpotsAndCentroids(dataTopDir, wellPath);
    results.experiment = repmat({'IAM942-34-subset'}, size(results, 1), 1);
    
    targetWellDir = [outExperimentDir, wName, '_scan'];
    
    if (~exist(targetWellDir, 'dir'))
        mkdir(targetWellDir)
    end
    
    writetable(results, [targetWellDir, '/centroidsAndTheirSpots.csv'])
    
    wellPath_stx = [experimentDir, wName, '_stacks'];
    if (~exist([outDir, wellPath_stx], 'dir'))
        mkdir([outDir, wellPath_stx])
    end
    extractStackSpotCountsAndMasks(dataTopDir, wellPath_stx, [outDir, wellPath_stx])
    
end

end