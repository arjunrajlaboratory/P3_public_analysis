function [] = IAMHCF30_subset_spotsAndCentroids(dataTopDir, outDir)


experimentDir = 'transdiff/201900629_IAMHCF30-subset/';

outExperimentDir = [outDir experimentDir];

if (~exist(outExperimentDir, 'dir'))
    mkdir(outExperimentDir)
end

readme = ReadYaml([dataTopDir, experimentDir, 'readMe.yaml']);

copyfile([dataTopDir, experimentDir, 'readMe.yaml'], [outDir, experimentDir]);

wellMap.wellName = strsplit(readme.wellMap.wellName, ',');
wellMap.condition = strsplit(readme.wellMap.condition, ',');


for i = 1:8
    
    wName = wellMap.wellName{i};
    wellPath = [experimentDir, wName, '_scan'];
    
    results = extractSpotsAndCentroids(dataTopDir, wellPath);
    results.experiment(:) = repmat({'IAMHCF30-subset'}, size(results, 1), 1);
    
    targetWellDir = [outExperimentDir, wName, '_scan'];
    
    if (~exist(targetWellDir, 'dir'))
        mkdir(targetWellDir)
    end
    
    writetable(results, [targetWellDir, '/centroidsAndTheirSpots.csv'])
    
%     wellPath_stx = [experimentDir, wName, '_stacks'];
%     if (~exist([outDir, wellPath_stx], 'dir'))
%         mkdir([outDir, wellPath_stx])
%     end
%     extractStackSpotCountsAndMasks(dataTopDir, wellPath_stx, [outDir, wellPath_stx])
    
end

end