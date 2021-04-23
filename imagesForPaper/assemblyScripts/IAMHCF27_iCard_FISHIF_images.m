function [outputArg1,outputArg2] = IAMHCF27_iCard_FISHIF_images(dataTopDir, outDir)
%UNTITLED Summary of this function goes here

%% contrast consistently
contrastAdjustValues.dapi = {[1200,10000]/65535};
contrastAdjustValues.alexa = {[684,2230]/65535};
contrastAdjustValues.cy = {[4000,12000]/65535};
contrastAdjustValues.tmr = {[883,2400]/65535};
contrastAdjustValues.gfp = {[2800,65300]/65535};


%% FISH-IF in iCard on 12-well punchout

stackNum = '002'; % pick the same tile for everything

% wellNum = 2;

R = [350 350 650 650]; % same part of each tile

slices = 1:7;

experimentDir = 'transdiff/20190528_HCF27_FISHIF/';

% projectPath = '/Volumes/IAM-BKUP4/cellid/transdiff/20190612_HCF-942-29-subset/942-29-subset/';

if ~exist([outDir, 'transdiff/20190528_HCF27_FISHIF/overlays/iCard'], 'dir')
    mkdir([outDir, 'transdiff/20190528_HCF27_FISHIF/overlays/iCard'])
end

copyfile([dataTopDir, experimentDir, 'readMe.yaml'], [outDir, 'transdiff/20190528_HCF27_FISHIF/']);

wellDir = 'iCard_stacks/';

fprintf(['working on overlays for ', wellDir, '\n']);

wellOut = [outDir, 'transdiff/20190528_HCF27_FISHIF/overlays/iCard/'];

cd([dataTopDir, 'transdiff/20190528_HCF27_FISHIF/iCard_stacks/'])

dapiStack = readmm(sprintf('dapi%s.tif', stackNum), slices);
tmrStack = readmm(sprintf('tmr%s.tif', stackNum), slices);
cyStack = readmm(sprintf('cy%s.tif', stackNum), slices);
alexaStack = readmm(sprintf('alexa%s.tif', stackNum), slices);
gfpStack = readmm(sprintf('gfp%s.tif', stackNum), slices);

images.channels =   {'dapi', 'alexa', 'cy', 'tmr', 'gfp'};
images.adj = {imadjust(imcrop(max(dapiStack.imagedata, [], 3), R), contrastAdjustValues.dapi{1}, []), ...
    imadjust(imcrop(max(alexaStack.imagedata, [], 3), R), contrastAdjustValues.alexa{1}, []), ...
    imadjust(imcrop(max(cyStack.imagedata, [], 3), R), contrastAdjustValues.cy{1}, []), ...
    imadjust(imcrop(max(tmrStack.imagedata, [], 3), R), contrastAdjustValues.tmr{1}, []),...
    imadjust(imcrop(max(gfpStack.imagedata, [], 3), R), contrastAdjustValues.gfp{1}, [])};

alexaDapiOverlay = cat(3, images.adj{2}, images.adj{2}, images.adj{2} + images.adj{1});
cyDapiOverlay = cat(3, images.adj{3}, images.adj{3}, images.adj{3} + images.adj{1});
tmrDapiOverlay = cat(3, images.adj{4}, images.adj{4}, images.adj{4} + images.adj{1});
gfpDapiOverlay = cat(3, images.adj{5}, images.adj{5}, images.adj{5} + images.adj{1});

imwrite(alexaDapiOverlay, fullfile(sprintf([wellOut 'alexa%s.tiff'], stackNum)));
imwrite(cyDapiOverlay, fullfile(sprintf([wellOut 'cy%s.tiff'], stackNum)));
imwrite(tmrDapiOverlay, fullfile(sprintf([wellOut 'tmr%s.tiff'], stackNum)));
imwrite(gfpDapiOverlay, fullfile(sprintf([wellOut 'gfp%s.tiff'], stackNum)));


end
