function [] = iCard942_repIF(dataTopDir, projectTopDir, imageSubdir)
%UNTITLED Summary of this function goes here

%% contrast consistently
contrastAdjustValues.dapi = {[1200,26000]/65535};
% contrastAdjustValues.alexa = {[684,2230]/65535};
% contrastAdjustValues.cy = {[4000,12000]/65535};
% contrastAdjustValues.tmr = {[883,2400]/65535};
contrastAdjustValues.gfp = {[4380,36000]/65535};


%% FISH-IF in iCard, just care about IF for now

stackNum = '003'; % pick the same tile for everything

% wellNum = 2;

R = [1 1 1023 1023]; % all of tile

slices = 7:20;

experimentDir = 'transdiff/20181116_testFISH/';

% projectPath = '/Volumes/IAM-BKUP4/cellid/transdiff/20190612_HCF-942-29-subset/942-29-subset/';
outFolder = [projectTopDir, imageSubdir, experimentDir];

if ~exist(outFolder, 'dir')
    mkdir(outFolder)
end

% copyfile([dataTopDir, experimentDir, 'readMe.yaml'], [outDir, 'transdiff/20190528_HCF27_FISHIF/']);

wellDir = 'GM942-iCard/w2/';

fprintf(['working on overlays for ', wellDir, '\n']);

wellOut = [outFolder, '/transdiff/20181116_testFISH/overlays/iCard/'];

cd([dataTopDir, 'transdiff/20181116_testFISH/GM942-iCard/w2/'])

dapiStack = readmm(sprintf('dapi%s.tif', stackNum), slices);
% tmrStack = readmm(sprintf('tmr%s.tif', stackNum), slices);
% cyStack = readmm(sprintf('cy%s.tif', stackNum), slices);
% alexaStack = readmm(sprintf('alexa%s.tif', stackNum), slices);
gfpStack = readmm(sprintf('gfp%s.tif', stackNum), slices);

images.channels =   {'dapi', 'gfp'};
images.adj = {imadjust(imcrop(max(dapiStack.imagedata, [], 3), R), contrastAdjustValues.dapi{1}, []), ...
%     imadjust(imcrop(max(alexaStack.imagedata, [], 3), R), contrastAdjustValues.alexa{1}, []), ...
%     imadjust(imcrop(max(cyStack.imagedata, [], 3), R), contrastAdjustValues.cy{1}, []), ...
%     imadjust(imcrop(max(tmrStack.imagedata, [], 3), R), contrastAdjustValues.tmr{1}, []),...
    imadjust(imcrop(max(gfpStack.imagedata, [], 3), R), contrastAdjustValues.gfp{1}, [])};

% alexaDapiOverlay = cat(3, images.adj{2}, images.adj{2}, images.adj{2} + images.adj{1});
% cyDapiOverlay = cat(3, images.adj{3}, images.adj{3}, images.adj{3} + images.adj{1});
% tmrDapiOverlay = cat(3, images.adj{4}, images.adj{4}, images.adj{4} + images.adj{1});
gfpDapiOverlay = cat(3, images.adj{2}, images.adj{2}, images.adj{2} + images.adj{1});

% imwrite(alexaDapiOverlay, fullfile(sprintf([wellOut 'alexa%s.tiff'], stackNum)));
% imwrite(cyDapiOverlay, fullfile(sprintf([wellOut 'cy%s.tiff'], stackNum)));
% imwrite(tmrDapiOverlay, fullfile(sprintf([wellOut 'tmr%s.tiff'], stackNum)));
imwrite(gfpDapiOverlay, fullfile(sprintf([outFolder 'gfp%s.tiff'], stackNum)));


end
