function [overlay] = extract4XWellOverlay(originalFile, ULcoord1, ULcoord2, contrastConds)
% grabs a 4X well-sized image at the indicated location in the indicated
% file
% well = wells{wid};

R = [ULcoord1 ULcoord2 384 384];

img = readmm(originalFile);

images.channels =   {'dapi', 'trans', 'tmr', 'gfp'};
images.adj = {imadjust(imcrop(max(img.imagedata, [], 3), R), contrastConds, [])};

overlay = cat(3, images.adj{1}, images.adj{1}, images.adj{1});

end

