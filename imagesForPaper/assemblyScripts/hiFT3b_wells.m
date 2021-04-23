function [overlay] = hiFT3b_wells(dataTopDir, projectTopDir, imageSubdir)

% extract representative images of 
% - scrambled shRNA control - plate 1 B4-6
% - CEBPB
% - PRRX2
% - RUNX1

% set contrast
contrastAdjustValues.trans = {[160,65000]/65535};
contrastAdjustValues.tmr = {[110,150]/65535};

experimentDir = 'hiFT/20181002_IAMhiFT3b/';

outFolder = [projectTopDir, imageSubdir, experimentDir];

if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

% scrambled shRNA 
transFile_SHC002 = [dataTopDir experimentDir 'Scan001_w1_s1015_t1.TIF'];
tmrFile_SHC002 = [dataTopDir experimentDir 'Scan001_w2_s1015_t1.TIF'];
UL1 = 1132;
UL2 = 398;

trans_SHC002 = extract4XWellOverlay(transFile_SHC002, UL1, UL2, contrastAdjustValues.trans{1});
tmr_SHC002 = extract4XWellOverlay(tmrFile_SHC002, UL1, UL2, contrastAdjustValues.tmr{1});

imwrite(trans_SHC002, [outFolder, 'SHC002_trans.tiff']);
imwrite(imcomplement(tmr_SHC002), [outFolder, 'SHC002_AlkPhos.tiff']);

% CEBPB
transFile_CEBPB = [dataTopDir experimentDir 'Scan001_w1_s1015_t1.TIF'];
tmrFile_CEBPB = [dataTopDir experimentDir 'Scan001_w2_s1015_t1.TIF'];
UL1 = 744;
UL2 = 1176;

trans_CEBPB = extract4XWellOverlay(transFile_CEBPB, UL1, UL2, contrastAdjustValues.trans{1});
tmr_CEBPB = extract4XWellOverlay(tmrFile_CEBPB, UL1, UL2, contrastAdjustValues.tmr{1});

imwrite(trans_CEBPB, [outFolder, 'CEBPB_trans.tiff']);
imwrite(imcomplement(tmr_CEBPB), [outFolder, 'CEBPB_AlkPhos.tiff']);

% PRRX2
transFile_PRRX2 = [dataTopDir experimentDir 'Scan002_w1_s1015_t1.TIF'];
tmrFile_PRRX2 = [dataTopDir experimentDir 'Scan002_w2_s1015_t1.TIF'];
UL1 = 1525;
UL2 = 742;

trans_PRRX2 = extract4XWellOverlay(transFile_PRRX2, UL1, UL2, contrastAdjustValues.trans{1});
tmr_PRRX2 = extract4XWellOverlay(tmrFile_PRRX2, UL1, UL2, contrastAdjustValues.tmr{1});

imwrite(trans_PRRX2, [outFolder, 'PRRX2_trans.tiff']);
imwrite(imcomplement(tmr_PRRX2), [outFolder, 'PRRX2_AlkPhos.tiff']);

% RUNX1
transFile_RUNX1 = [dataTopDir experimentDir 'Scan002_w1_s1015_t1.TIF'];
tmrFile_RUNX1 = [dataTopDir experimentDir 'Scan002_w2_s1015_t1.TIF'];
UL1 = 1905;
UL2 = 1127;

trans_RUNX1 = extract4XWellOverlay(transFile_RUNX1, UL1, UL2, contrastAdjustValues.trans{1});
tmr_RUNX1 = extract4XWellOverlay(tmrFile_RUNX1, UL1, UL2, contrastAdjustValues.tmr{1});

imwrite(trans_RUNX1, [outFolder, 'RUNX1_trans.tiff']);
imwrite(imcomplement(tmr_RUNX1), [outFolder, 'RUNX1_AlkPhos.tiff']);
end