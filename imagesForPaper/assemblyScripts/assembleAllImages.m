%% assembleAllImages.m - create all image panels for paper
% assembles all images for panels in the paper
% transdifferentiation images
% - smFISH scan centroid maps for all replicates
% - representative smFISH image overlays from stacks
% - smFISH/cTnT IF positive controls in iPSC-derived cardiomyocytes
% includes hiF-T
% - Tra-1-60/AlkPhos overlays
% - representative AlkPhos-stained wells

%% Data and project folders
% adjust these based on where you have stored the data and where you want
% the output to be stored

% top-level location of all data to be processed into image panels
% dataTopDir = '/Volumes/IAM-BKUP4/cellid/';
% tempDataTopDir = '/Volumes/IAM-BKUP1/cellid/';
dataTopDir = '~/Dropbox (RajLab)/Shared_IanM/cellid_201807_onward/rawData/';

% top-level location in which image panels will be stored
projectTopDir = '~/Dropbox (RajLab)/Shared_IanM/cellid_201807_onward/';

% subdirectory of projectTopDir in which to place processed images
imageSubDir = 'processedImages_test430/';

if ~exist([projectTopDir, imageSubDir], 'dir')
    mkdir([projectTopDir, imageSubDir])
end

%% example images and videos from perturbation experiments
% fibroblasts
RNAtagTest5_repImg(dataTopDir, [projectTopDir, imageSubDir])

% iCards
iCardTest5_repImgVideo(dataTopDir, [projectTopDir, imageSubDir])

%% cardiac transdifferentiation experiments - smFISH and cTnT IF
% for each experiment, extract centroid locations and spot counts in
% representative wells

% GM00942 dermal fibroblasts - replicate 1
IAM942_34_subset_spotsAndCentroids(dataTopDir, [projectTopDir, imageSubDir])
IAM942_34_subset_images(dataTopDir, [projectTopDir, imageSubDir])
% IAM942_34_subset_stacksImages(dataTopDir, [projectTopDir, imageSubDir])

% GM00942 dermal fibroblasts - replicate 2
IAM942_30_subset_spotsAndCentroids(dataTopDir, [projectTopDir, imageSubDir])

% GM11169 cardiac fibroblasts - replicate 1
IAM11169_30_subset_spotsAndCentroids(dataTopDir, [projectTopDir, imageSubDir])

% GM11169 cardiac fibroblasts - replicate 2
IAM11169_28_subset_spotsAndCentroids(dataTopDir, [projectTopDir, imageSubDir])

% immortalized cardiac fibroblast line - replicate 1
IAMHCF30_subset_spotsAndCentroids(dataTopDir, [projectTopDir, imageSubDir])

% immortalized cardiac fibroblast line - replicate 2
IAMHCF22_subset_spotsAndCentroids(dataTopDir, [projectTopDir, imageSubDir])

% iPSC-derived cardiomyocyte FISH-IF controls
IAMHCF27_iCard_FISHIF_images(dataTopDir, [projectTopDir, imageSubDir])

% representative iPSC-derived cardiomyocyte IF


%% hiF-T-iPSC Alk Phos calorimetry and Tra-1-60 IF
% stores images in experiment-specific subdirectory of 
% [projectTopDir imageSubDir] called 'hiFT/'

% Tra-1-60 IF / Alk Phos fluorescent colorimetry correlation
Tra_1_60_AlkPhos_overlays(dataTopDir, projectTopDir, imageSubDir)

% Representative Alk Phos fluorescent colorimetry for colony counts
hiFT4b_wells(dataTopDir, projectTopDir, imageSubDir)
hiFT3b_wells(dataTopDir, projectTopDir, imageSubDir)

