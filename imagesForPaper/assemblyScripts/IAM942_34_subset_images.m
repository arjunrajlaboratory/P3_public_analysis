function [outputArg1,outputArg2] = IAM942_34_subset_images(dataTopDir, outDir)
%UNTITLED Summary of this function goes here

%% contrast consistently
contrastAdjustValues.dapi = {[1200,10000]/65535};
contrastAdjustValues.alexa = {[625,1750]/65535};
contrastAdjustValues.cy = {[3180,6400]/65535};
contrastAdjustValues.tmr = {[800,1900]/65535};
contrastAdjustValues.gfp = {[1400,2500]/65535};


%% TNNT2+ cell? in 6UP-only well 5 of 942-29-subset
% Script to create composite zstack of NGFR-mNG2 tagged WM989 A6-G3 C10-C2 clone E9 RNA FISH images. 
% Uses readmm.m function from Rajlabimagetools. Script can be run from the
% project path. 

stackNum = '005'; % pick the same tile for everything

wellNum = 2;

R = [268 24 600 400]; % same part of each tile

slices = 3:36;

% projectPath = '/Volumes/IAM-BKUP4/cellid/transdiff/20190612_HCF-942-29-subset/942-29-subset/';

if ~exist([outDir, 'transdiff/20190731_IAM942-34-subset/overlays'], 'dir')
    mkdir([outDir, 'transdiff/20190731_IAM942-34-subset/overlays'])
end

wellDir = ['well', num2str(wellNum), '_stacks/'];

fprintf(['working on overlays for ', wellDir, '\n']);

wellOut = [outDir, 'transdiff/20190731_IAM942-34-subset/overlays/'];

cd([dataTopDir, 'transdiff/20190731_IAM942-34-subset/well', num2str(wellNum), '_stacks/'])

dapiStack = readmm(sprintf('dapi%s.tif', stackNum), slices);
tmrStack = readmm(sprintf('tmr%s.tif', stackNum), slices);
cyStack = readmm(sprintf('cy%s.tif', stackNum), slices);
alexaStack = readmm(sprintf('alexa%s.tif', stackNum), slices);

images.channels =   {'dapi', 'alexa', 'cy', 'tmr'};
images.adj = {imadjust(imcrop(max(dapiStack.imagedata, [], 3), R), contrastAdjustValues.dapi{1}, []), ...
    imadjust(imcrop(max(alexaStack.imagedata, [], 3), R), contrastAdjustValues.alexa{1}, []), ...
    imadjust(imcrop(max(cyStack.imagedata, [], 3), R), contrastAdjustValues.cy{1}, []), ...
    imadjust(imcrop(max(tmrStack.imagedata, [], 3), R), contrastAdjustValues.tmr{1}, [])};

alexaDapiOverlay = cat(3, images.adj{2}, images.adj{2}, images.adj{2} + images.adj{1});
cyDapiOverlay = cat(3, images.adj{3}, images.adj{3}, images.adj{3} + images.adj{1});
tmrDapiOverlay = cat(3, images.adj{4}, images.adj{4}, images.adj{4} + images.adj{1});

imwrite(alexaDapiOverlay, fullfile(sprintf([wellOut 'alexa%s.tiff'], stackNum)));
imwrite(cyDapiOverlay, fullfile(sprintf([wellOut 'cy%s.tiff'], stackNum)));
imwrite(tmrDapiOverlay, fullfile(sprintf([wellOut 'tmr%s.tiff'], stackNum)));


end

