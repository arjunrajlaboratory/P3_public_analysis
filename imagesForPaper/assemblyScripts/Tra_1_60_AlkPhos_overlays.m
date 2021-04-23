function [] = Tra_1_60_AlkPhos_overlays(dataTopDir, projectTopDir, imageSubdir)

%% contrast from scans in well 2 (untransduced hiFT-iPSC)
contrastAdjustValues.dapi = {[1200,10000]/65535};
contrastAdjustValues.alexa = {[700,1100]/65535};
contrastAdjustValues.cy = {[2400,4600]/65535};
contrastAdjustValues.tmr = {[670,1300]/65535};
contrastAdjustValues.gfp = {[800,2000]/65535};
contrastAdjustValues.trans = {[1200,65000]/65535};

%%

cd([dataTopDir, 'hiFT/20190727-hiFT3a_Tra-1-60/'])

wells = {'well1_PO_noTra', 'well2_PO_withTra', 'well3_SHC002_withTra',...
    'well4_CEBPBsh_withTra', 'well5_PRRX2sh_withTra', 'well6_RUNX1sh_withTra'};

colonyUL = {{10650, 2680},...
    {5650, 3020},...
    {7220, 1340},...
    {3080, 6450},...
    {6480, 10260},...
    {3040, 3540}};

outFolder = [projectTopDir, imageSubdir, 'hiFT/Tra-1-60_AlkPhos/'];

if ~exist([projectTopDir, imageSubdir, 'hiFT/'], 'dir')
    mkdir([projectTopDir, imageSubdir, 'hiFT/'])
end

if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

for wid = 1:6
    
    well = wells{wid};
    
    UL_x = colonyUL{wid}{1} - 500;
    UL_y = colonyUL{wid}{2} - 500;
    
    R = [UL_x UL_y 1000 1000];
    
    bigimg = imread([well, '.tif']);
    
    transSmall = imresize(bigimg(:,:,1), 0.1);
    dapiSmall = imresize(bigimg(:,:,2), 0.1);
    tmrSmall = imresize(bigimg(:,:,3), 0.1);
    gfpSmall = imresize(bigimg(:,:,4), 0.1);
    
    images.channels =   {'dapi', 'trans', 'tmr', 'gfp'};
    images.adj = {imadjust(dapiSmall, contrastAdjustValues.dapi{1}, []), ...
        imadjust(transSmall, contrastAdjustValues.trans{1}, []), ...
        imadjust(tmrSmall, contrastAdjustValues.tmr{1}, []), ...
        imadjust(gfpSmall, contrastAdjustValues.gfp{1}, []), ...
        imadjust(imcrop(bigimg(:,:,1), R), contrastAdjustValues.trans{1}, []),...
        imadjust(imcrop(bigimg(:,:,3), R), contrastAdjustValues.tmr{1}, []),...
        imadjust(imcrop(bigimg(:,:,4), R), contrastAdjustValues.gfp{1}, [])};
    
    transDapiOverlay = cat(3, images.adj{2}, images.adj{2}, images.adj{2} + images.adj{1});
    tmrDapiOverlay = cat(3, images.adj{3}, images.adj{3}, images.adj{3} + images.adj{1});
    gfpDapiOverlay = cat(3, images.adj{4}, images.adj{4}, images.adj{4} + images.adj{1});
    
    imwrite(transDapiOverlay, fullfile(sprintf([outFolder '%s_trans_small.tiff'], well)));
    imwrite(tmrDapiOverlay, fullfile(sprintf([outFolder '%s_tmr_small.tiff'], well)));
    imwrite(gfpDapiOverlay, fullfile(sprintf([outFolder '%s_gfp_small.tiff'], well)));
    
    transAlone = cat(3, images.adj{2}, images.adj{2}, images.adj{2});
    tmrRed = cat(3, images.adj{3}, zeros(size(images.adj{3})), zeros(size(images.adj{3})));
    gfpGreen = cat(3, zeros(size(images.adj{4})), images.adj{4}, zeros(size(images.adj{4})));
%     tmrgfpOverlay = cat(3, images.adj{3}, images.adj{4}, zeros(size(images.adj{4})));
    colonyAlone = cat(3, images.adj{5}, images.adj{5}, images.adj{5});
    colonyTMR = cat(3, images.adj{6}, zeros(size(images.adj{6})), zeros(size(images.adj{6})));
    colonyGFP = cat(3, zeros(size(images.adj{7})), images.adj{7}, zeros(size(images.adj{7})));
    
    imwrite(transAlone, fullfile(sprintf([outFolder '%s_trans_alone_small.tiff'], well)));
    imwrite(tmrRed, fullfile(sprintf([outFolder '%s_tmr_red_small.tiff'], well)));
    imwrite(gfpGreen, fullfile(sprintf([outFolder '%s_gfp_green_small.tiff'], well)));
    imwrite(colonyAlone, fullfile(sprintf([outFolder '%s_colony_alone_small.tiff'], well)));
    imwrite(colonyTMR, fullfile(sprintf([outFolder '%s_colony_TMR_small.tiff'], well)));
    imwrite(colonyGFP, fullfile(sprintf([outFolder '%s_colony_GFP_small.tiff'], well)));
    
end
clear bigimg

end

