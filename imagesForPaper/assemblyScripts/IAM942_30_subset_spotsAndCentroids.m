function [] = IAM942_30_subset_spotsAndCentroids(dataTopDir, outDir)


experimentDir = 'transdiff/20190701_IAM11169-942-30-subset/942-30-subset/';
experimentOutDir = 'transdiff/20190701_IAM942-30-subset/';
outExperimentDir = [outDir experimentOutDir];

if (~exist(outExperimentDir, 'dir'))
    mkdir(outExperimentDir)
end

readme = ReadYaml([dataTopDir, experimentDir, 'readMe.yaml']);

copyfile([dataTopDir, experimentDir, 'readMe.yaml'], [outDir, experimentOutDir]);

wellMap.wellName = strsplit(readme.wellMap.wellName, ',');
wellMap.condition = strsplit(readme.wellMap.condition, ',');


for i = 1:6
    
    wName = wellMap.wellName{i};
    wellPath = [experimentDir, wName, '_scan'];
    
    results = extractSpotsAndCentroids(dataTopDir, wellPath);
%     results.experiment(:) = repmat({'IAM942-30-subset'}, size(results,
%     1), 1); % run this line instead of 28 if in MATLAB_2019a
    results.experiment = repmat({'IAM942-30-subset'}, size(results, 1), 1);
    
    targetWellDir = [outExperimentDir, wName, '_scan'];
    
    if (~exist(targetWellDir, 'dir'))
        mkdir(targetWellDir)
    end
    
    writetable(results, [targetWellDir, '/centroidsAndTheirSpots.csv'])
%     
%     wellPath_stx = [experimentDir, wName, '_stacks'];
%     if (~exist([outDir, wellPath_stx], 'dir'))
%         mkdir([outDir, wellPath_stx])
%     end
%     extractStackSpotCountsAndMasks(dataTopDir, wellPath_stx, [outDir, wellPath_stx])
    
end

end