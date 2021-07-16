clear all; close all; clc;
%% copy images 
outdir = 'allImages';
if ~exist(outdir)
    mkdir(outdir);
end

conds = {'ESI','LKO'}; %conditions, subfolder name in a master data folder
times = {'20','30'};
ncols = 6;
%%
q = 1; 
for ii = 1:length(conds)
    for jj = 1:length(times)
        for kk = 1:ncols
            filename = [conds{ii} times{jj} 'h' filesep 'colony' int2str(kk) '.tif'];
            outfilename = [outdir filesep 'colony' sprintf('%.2d',q) '.tif'];
            copyfile(filename,outfilename);
            q = q + 1;
        end
        allConds{(q-1)/ncols} = [conds{ii} times{jj} 'h'];
        filenrs{(q-1)/ncols} =  (q-6):(q-1);
    end
end
%% 

%%
addpath(genpath('/Documents/ImageProcessing-master')); 
warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');

mainDataDir = './';
dataDir =  fullfile('allImages');
mdatafile = '//Volumes/Li/20200826b/allData/metadata.mat';
load(mdatafile);
fInFormat = 'colony%.2d.tif';
filename = sprintf(fInFormat, 1);
colRadius = 310;

%meta = MetadataMicropattern(fullfile(dataDir,filename));

% fix some stuff manually
 meta.xres = 0.625;
 meta.yres = 0.625;
% 
 meta.colRadiiMicron = colRadius;
 meta.colRadiiPixel = round(meta.colRadiiMicron/meta.xres);
% 
% fileRange = folderFilesFromKeyword(dataDir,'Image');

%BLOCK6 = {{1 5} {11 14} {19 22} {28 32}};  %Block array / removed a bad photo N0. 15
meta.channelLabel = {'DAPI','LeftyRNA','BraRNA','Axin2RNA'};

%filenrs = {25:30};
% meta.channelLabel = {'DAPI','SIX1','PAX3','SNAI1'};

meta.conditions = allConds;

save(fullfile(dataDir,'metaData.mat'),'meta');

DAPIChannel = 1;

%% process (save DAPI channel MIP for segmentation)

% although findColonies can find multiple colonies in a single image,
% for the LSM there is one per image, which this code below assumes
findColParam = struct('sclose', 6, 'sopen', 8, 'checkcontained', false,...
                            'minArea', [],'convhull', true);
close all;
colonies(numel([filenrs{:}])) = Colony;


for coli = [filenrs{:}]
  
    % cleanScale in micron
    param = {'DAPIChannel',DAPIChannel, 'colID',coli, 'adjustmentFactor', 0.2,'clparameters',findColParam,'thresh',[-1 -1 -1 -1]};
    %  thresh param, set to the maximum value you want to allow. 1 number per channel. set to -1 not to apply to that channel.
    filename = sprintf(fInFormat, coli);
    colony = processOneColonyImage(filename, dataDir, param);
    colonies(coli) = colony;
    colonies(coli).setID(coli);
    % setID overwrites filename, which was the right thing for epi, not
    % here
    colonies(coli).filename = filename; 
end

save(fullfile(dataDir,'colonies'), 'colonies');

%% show plot of different conditions side by side

load(fullfile(dataDir,'colonies'), 'colonies');

% first normalize by DAPI, then scale all profiles from 0 to 1
doubleNormalize = true; 

for i = 1:numel(meta.conditions)   
%for i = 1:5
    coloniesCombined{i} = colonies(filenrs{i});
end

figure('Position',[0 0 900 900]);
plotMultipleAveragesNoSegmentation(meta, colRadius, DAPIChannel,...
                                coloniesCombined, meta.conditions, doubleNormalize)
saveas(gcf,fullfile(dataDir,'radialProfilesCombined.png'));                        


%% show plot of different conditions side by side w/o DAPI  normalize

doubleNormalize = true;
figure('Position',[0 0 900 900]);
plotMultipleAveragesNoSegmentation(meta, colRadius, [],...
                                coloniesCombined, meta.conditions, doubleNormalize)
saveas(gcf,fullfile(dataDir,'radialProfilesCombinedNoDAPI.png'));                        

%% overlay conditions per stain

figure('Position',[0 0 2000 2000]);
plotConditionOverlayNoSegmentation(meta, colRadius, [],...
                                coloniesCombined(8:14), meta.conditions(8:14), doubleNormalize);
saveas(gcf,fullfile(dataDir,'radialProfilesConditionOverlay.png'));
%%
cleanupOpts = struct('minArea', 30, 'cytoplasmicLevels', false);
opts = struct('cleanupOptions', cleanupOpts);
colDir = dataDir; 

for coli = [filenrs{:}]
    
    % extract segmented data
    colonies(coli).extractData(colDir, DAPIChannel, opts);

    % radial binning of segmented data
    colType = find(meta.colRadiiMicron == colonies(coli).radiusMicron);
    colonies(coli).makeRadialAvgSeg()
end

save(fullfile(dataDir,'colonies'), 'colonies');

%% load colonies and show side by side

load(fullfile(dataDir,'metaData.mat'),'meta');
load(fullfile(dataDir,'colonies'), 'colonies');
colDir = dataDir; 
% representative colonies for each condition
repcols = [1 2 3 4 5 6; 7 8 9 10 11 12; 13 14 15 16 17 18]'; 

ci = 1; % sox9

figure('Position',[0 0 1600 1000]),
m = 6;
n = numel(filenrs);

% set LUT from image of my choice, for SOX9 the 24h iwp is the brightest
col = colonies(filenrs{end}(1));
img = col.loadImage(colDir, ci);
img = max(img,[],3);
Ilim = round(stretchlim(img)*(2^16-1));

for i = 1:n
    for j = 1:m
   
        ii = (j-1)*n + i;

        col = colonies(filenrs{i}(m));
        colDir = dataDir;
        img = col.loadImage(colDir, ci);
        img = max(img,[],3);
        b= col.boundingBox;
        img = img(b(3):b(4),b(1):b(2));
        %bg = imopen(img, strel('disk',15));
        %img = img - bg;
%         if j==1 && i == 1
%             Ilim = round(stretchlim(img)*(2^16-1));
%         end
        subplot_tight(m,n,ii)
        imshow(img,Ilim);
        title(meta.conditions{i})
    end
end
saveas(gcf,fullfile(dataDir,['compareColonies_' meta.channelLabel{ci} '.png']));

