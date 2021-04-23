function [overlay1, overlay2, video1 ] = iCardTest5_repImgVideo(dataTopDir, projectTopDir, imageSubdir)

% extract representative images and video of 
% - P3 iCards after drugging

% set contrast
contrastAdjustValues.trans = {[12000,36000]/65535};
contrastAdjustValues.diff = {[0,2000]/65535};

experimentDir = 'iCard-942/videos/20170827_iCardTest5set100x_94h/';

outFolder = [projectTopDir, imageSubdir, experimentDir];

if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

% iCardTest5 well F7, WHI-P154 - drugged 
transFile = [dataTopDir experimentDir 'F7_001.TIF'];
UL1 = 260;
UL2 = 130;

R = [UL1 UL2 384 384];

img = readmm(transFile);

images1.channels =   {'trans'};
images1.adj = {imadjust(imcrop(img.imagedata(:,:,29), R), contrastAdjustValues.trans{:}, [])};
images2.channels =   {'trans'};
images2.adj = {imadjust(imcrop(img.imagedata(:,:,30), R), contrastAdjustValues.trans{:}, [])};

imdiff.channels =   {'trans'};
imdiff.adj = {imadjust(abs(imcrop(img.imagedata(:,:,30), R) - imcrop(img.imagedata(:,:,29), R)), contrastAdjustValues.diff{:}, [])};

% subset the video
i = 1;
video1.adj{1}(:,:,i) = imadjust(imcrop(img.imagedata(:,:,i), R), contrastAdjustValues.trans{:}, []);
imwrite(video1.adj{1}(:,:,i), [outFolder, 'F7_video_perturb.tiff']);
for i = 2:41
    
    video1.adj{1}(:,:,i) = imadjust(imcrop(img.imagedata(:,:,i), R), contrastAdjustValues.trans{:}, []);
    imwrite(video1.adj{1}(:,:,i), [outFolder, 'F7_video_perturb.tiff'], 'WriteMode','append');
end

% finish the stills
overlay1 = cat(3, images1.adj{1}, images1.adj{1}, images1.adj{1});
overlay2 = cat(3, images2.adj{1}, images2.adj{1}, images2.adj{1});

imwrite(overlay1, [outFolder, 'F7_still1_perturb.tiff']);
imwrite(overlay2, [outFolder, 'F7_still2_perturb.tiff']);

overlay3 = cat(3, images1.adj{1}, images2.adj{1}, zeros(size(images1.adj{1})));

imwrite(overlay3, [outFolder, 'F7_redgreen_perturb.tiff']);

overlay4 = cat(3, imdiff.adj{1}, imdiff.adj{1}, imdiff.adj{1});

imwrite(overlay4, [outFolder, 'F7_diff_perturb.tiff']);


% iCardTest5 well G9, DMSOonly - control 
transFile = [dataTopDir experimentDir 'G9_001.TIF'];
UL1 = 260;
UL2 = 130;

R = [UL1 UL2 384 384];

img = readmm(transFile);

images1.channels =   {'trans'};
images1.adj = {imadjust(imcrop(img.imagedata(:,:,29), R), contrastAdjustValues.trans{:}, [])};
images2.channels =   {'trans'};
images2.adj = {imadjust(imcrop(img.imagedata(:,:,30), R), contrastAdjustValues.trans{:}, [])};

imdiff.channels =   {'trans'};
imdiff.adj = {imadjust(abs(imcrop(img.imagedata(:,:,30), R) - imcrop(img.imagedata(:,:,29), R)), contrastAdjustValues.diff{:}, [])};

% subset the video
i = 1;
video1.adj{1}(:,:,i) = imadjust(imcrop(img.imagedata(:,:,i), R), contrastAdjustValues.trans{:}, []);
imwrite(video1.adj{1}(:,:,i), [outFolder, 'G9_video_control.tiff']);
for i = 2:41
    
    video1.adj{1}(:,:,i) = imadjust(imcrop(img.imagedata(:,:,i), R), contrastAdjustValues.trans{:}, []);
    imwrite(video1.adj{1}(:,:,i), [outFolder, 'G9_video_control.tiff'], 'WriteMode','append');
end

% finish the stills
overlay1 = cat(3, images1.adj{1}, images1.adj{1}, images1.adj{1});
overlay2 = cat(3, images2.adj{1}, images2.adj{1}, images2.adj{1});

imwrite(overlay1, [outFolder, 'G9_still1_control.tiff']);
imwrite(overlay2, [outFolder, 'G9_still2_control.tiff']);

overlay3 = cat(3, images1.adj{1}, images2.adj{1}, zeros(size(images1.adj{1})));

imwrite(overlay3, [outFolder, 'G9_redgreen_control.tiff']);

overlay4 = cat(3, imdiff.adj{1}, imdiff.adj{1}, imdiff.adj{1});

imwrite(overlay4, [outFolder, 'G9_diff_control.tiff']);


end