function [overlay] = RNAtagTest5_repImg(dataTopDir, projectTopDir, imageSubdir)

% extract representative images of 
% - P3 fibroblasts after drugging
%%
% set contrast
contrastAdjustValues.trans = {[3600,14000]/65535};

experimentDir = 'GM00942/scans/20170817_RNAtagTest5_drugged_96h_tallScan/';

outFolder = [projectTopDir, imageSubdir, experimentDir];

if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

% RNAtagTest5 well 15, WHI-P154 - drugged 
transFile = [dataTopDir experimentDir 'trans015.tif'];
UL1 = 25;
UL2 = 476;

R = [UL1 UL2 384 384];

img = readmm(transFile);

images.channels =   {'trans'};
images.adj = {imadjust(imcrop(img.imagedata(:,:,15), R), contrastAdjustValues.trans{:}, [])};

overlay = cat(3, images.adj{1}, images.adj{1}, images.adj{1});
imwrite(overlay, [outFolder, 'trans015_perturb.tiff']);

% RNAtagTest5 well 25, DMSO-only control 
transFile = [dataTopDir experimentDir 'trans025.TIF'];
UL1 = 25;
UL2 = 476;

R = [UL1 UL2 384 384];

img = readmm(transFile);

images.channels =   {'trans'};
images.adj = {imadjust(imcrop(img.imagedata(:,:,15), R), contrastAdjustValues.trans{:}, [])};

overlay = cat(3, images.adj{1}, images.adj{1}, images.adj{1});
imwrite(overlay, [outFolder, 'trans025_control.tiff']);
end