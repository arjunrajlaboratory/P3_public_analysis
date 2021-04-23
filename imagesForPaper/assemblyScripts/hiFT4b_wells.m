function [overlay] = hiFT4b_wells_forPoster(dataTopDir, projectTopDir, imageSubdir)

% extract representative images of 
% - scrambled shRNA control - plate 1 B4-6
% - ZNF652 - plate 5 C4-6
% - ATOH8 - plate 5 A4-6
% - KLF13 - plate 3 B1-6
% - ZBTB38 - plate 4 D1-3

% set contrast
contrastAdjustValues.trans = {[160,65000]/65535};
% contrastAdjustValues.tmr = {[120,162]/65535};
contrastAdjustValues.tmr = {[120,180]/65535};

experimentDir = 'hiFT/20181108_IAMhiFT4b/';

outFolder = [projectTopDir, imageSubdir, experimentDir];

if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

% scrambled shRNA 
transFile_SHC002 = [dataTopDir experimentDir 'Scan001_w1_s1081_t1.TIF'];
tmrFile_SHC002 = [dataTopDir experimentDir 'Scan001_w2_s1081_t1.TIF'];
UL1 = 1124;
UL2 = 408;

trans_SHC002 = extract4XWellOverlay(transFile_SHC002, UL1, UL2, contrastAdjustValues.trans{1});
tmr_SHC002 = extract4XWellOverlay(tmrFile_SHC002, UL1, UL2, contrastAdjustValues.tmr{1});

imwrite(trans_SHC002, [outFolder, 'SHC002_trans.tiff']);
imwrite(imcomplement(tmr_SHC002), [outFolder, 'SHC002_AlkPhos.tiff']);

% ZNF652
transFile_ZNF652 = [dataTopDir experimentDir 'Scan005_w1_s1081_t1.TIF'];
tmrFile_ZNF652 = [dataTopDir experimentDir 'Scan005_w2_s1081_t1.TIF'];
UL1 = 1124;
UL2 = 778;

trans_ZNF652 = extract4XWellOverlay(transFile_ZNF652, UL1, UL2, contrastAdjustValues.trans{1});
tmr_ZNF652 = extract4XWellOverlay(tmrFile_ZNF652, UL1, UL2, contrastAdjustValues.tmr{1});

imwrite(trans_ZNF652, [outFolder, 'ZNF652_trans.tiff']);
imwrite(imcomplement(tmr_ZNF652), [outFolder, 'ZNF652_AlkPhos.tiff']);

% ATOH8
transFile_ATOH8 = [dataTopDir experimentDir 'Scan005_w1_s1081_t1.TIF'];
tmrFile_ATOH8 = [dataTopDir experimentDir 'Scan005_w2_s1081_t1.TIF'];
UL1 = 1124;
UL2 = 42;

trans_ATOH8 = extract4XWellOverlay(transFile_ATOH8, UL1, UL2, contrastAdjustValues.trans{1});
tmr_ATOH8 = extract4XWellOverlay(tmrFile_ATOH8, UL1, UL2, contrastAdjustValues.tmr{1});

imwrite(trans_ATOH8, [outFolder, 'ATOH8_trans.tiff']);
imwrite(imcomplement(tmr_ATOH8), [outFolder, 'ATOH8_AlkPhos.tiff']);

% KLF13
transFile_KLF13 = [dataTopDir experimentDir 'Scan003_w1_s1081_t1.TIF'];
tmrFile_KLF13 = [dataTopDir experimentDir 'Scan003_w2_s1081_t1.TIF'];
UL1 = 755;
UL2 = 405;

trans_KLF13 = extract4XWellOverlay(transFile_KLF13, UL1, UL2, contrastAdjustValues.trans{1});
tmr_KLF13 = extract4XWellOverlay(tmrFile_KLF13, UL1, UL2, contrastAdjustValues.tmr{1});

imwrite(trans_KLF13, [outFolder, 'KLF13_trans.tiff']);
imwrite(imcomplement(tmr_KLF13), [outFolder, 'KLF13_AlkPhos.tiff']);

% ZBTB38
transFile_ZBTB38 = [dataTopDir experimentDir 'Scan004_w1_s1081_t1.TIF'];
tmrFile_ZBTB38 = [dataTopDir experimentDir 'Scan004_w2_s1081_t1.TIF'];
UL1 = 744;
UL2 = 1158;

trans_ZBTB38 = extract4XWellOverlay(transFile_ZBTB38, UL1, UL2, contrastAdjustValues.trans{1});
tmr_ZBTB38 = extract4XWellOverlay(tmrFile_ZBTB38, UL1, UL2, contrastAdjustValues.tmr{1});

imwrite(trans_ZBTB38, [outFolder, 'ZBTB38_trans.tiff']);
imwrite(imcomplement(tmr_ZBTB38), [outFolder, 'ZBTB38_AlkPhos.tiff']);

% CERS2
transFile_CERS2 = [dataTopDir experimentDir 'Scan003_w1_s1081_t1.TIF'];
tmrFile_CERS2 = [dataTopDir experimentDir 'Scan003_w2_s1081_t1.TIF'];
UL1 = 370;
UL2 = 1150;

trans_CERS2 = extract4XWellOverlay(transFile_CERS2, UL1, UL2, contrastAdjustValues.trans{1});
tmr_CERS2 = extract4XWellOverlay(tmrFile_CERS2, UL1, UL2, contrastAdjustValues.tmr{1});

imwrite(trans_CERS2, [outFolder, 'CERS2_trans.tiff']);
imwrite(imcomplement(tmr_CERS2), [outFolder, 'CERS2_AlkPhos.tiff']);


end