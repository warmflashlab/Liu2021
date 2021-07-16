
%%
%% -------------------- 1) specify filepaths, wells, channels details----------------------------
close all; clearvars; clc;
addpath(genpath('/Users/ll/Desktop/MatlabFunctions'));
%%
masterFolder = '/Volumes/WD Element/20201003/allData';
nWells = 6;
wellNames = {'ESI6h','LKO6h', 'ESI24h', 'LKO24h', 'ESI48hh', 'LKO48h',...
     'ESI','L1KOA','L2KOC'};
nColoniesPerWell = [6 6 6 6 6 6 6 6 6];% number of colonies per well
nImagesPerColonyX = 2; nImagesPerColonyY = 2; % number of images per colony in the X and Y direction
%%
rawDataFile1 = '20201003_f0000_w0000_z0000.tif';
rawDataFolder = [masterFolder filesep 'rawData'];
meta = MetadataMicropattern([rawDataFolder filesep rawDataFile1]);
%%
meta.colRadiiMicron = 350;
meta.channelNames = {'DAPI', 'Axin2','smad2',  'Bra' };
meta.channelLabel = {'nuclear', 'nuclear','nuclear',  'nuclear'};%[spots for RNA, nonMembrane for betacatenin]
meta.nChannels = numel(meta.channelNames);
save([masterFolder filesep 'metadata.mat'], 'meta');
%%
% overalap in pixels for stitching
overlapLR = 600;
overlapTB = 700;

%%
%% ---------------------- 2) read images, make maxZ projections, backgroundSubtract, --------------
%% ---------------------------------align images, make colonyMasks, save --------------------------


rawDataFiles = readAndorDirectory(rawDataFolder);
rawDataFiles = rawDataFiles(1);

processedDataFolder = [masterFolder filesep 'processedData'];
mkdir(processedDataFolder);
%%

colonyCounter = 1; % absolute colony counter

dapiChannel = find(cell2mat(cellfun(@(c)strcmp(c,'DAPI'),upper(meta.channelNames),'UniformOutput',false)));
membraneChannel = find(cell2mat(cellfun(@(c)strcmp(c,'NONMEMBRANE'),upper(meta.channelLabel),'UniformOutput',false)));
spotChannels = find(cell2mat(cellfun(@(c)strcmp(c,'SPOTS'),upper(meta.channelLabel),'UniformOutput',false)));
bfChannel = find(cell2mat(cellfun(@(c)strcmp(c,'BF'),upper(meta.channelNames),'UniformOutput',false)));


maxZChannels = setxor(rawDataFiles.w+1, bfChannel)-1;

nImagesPerColony = nImagesPerColonyX*nImagesPerColonyY;

for ii = 1:nWells
    tic;
    ii
    close all;
    wellFolder = [processedDataFolder filesep wellNames{ii}];
    mkdir(wellFolder);
    
    for jj = 1:nColoniesPerWell(ii)
        
        % make max z projections
        % background subtract
        imageIds = [(colonyCounter-1)*nImagesPerColony+1:(colonyCounter-1)*nImagesPerColony+nImagesPerColony]-1;
        images = cell(1,nImagesPerColony);
        imgCounter = 1;
        for kk = imageIds
            channelCounter = 1;
            for ll = maxZChannels
                maxImg = andorMaxIntensity(rawDataFiles,kk,0,ll); %[files,pos,time,chan]
                images{imgCounter}(:,:,channelCounter) = smoothAndBackgroundSubtractOneImage(maxImg);
                channelCounter = channelCounter+1;
            end
            % add the first z slice of brightfield channel as the last image of images array
            if ~isempty(bfChannel)
                imageName = getAndorFileName(rawDataFiles,kk,0,0,bfChannel-1);
                images{imgCounter}(:,:,channelCounter) = imread(imageName);
            end
            
            imgCounter = imgCounter+1;
        end
        
        if nImagesPerColony > 1
            finalImage = makeColonyImage(images, nImagesPerColonyX, nImagesPerColonyY, overlapLR, overlapTB); % align
            smoothMask = makeColonyMaskUsingDapiImage(finalImage(:,:,dapiChannel), meta); % make colony mask
            % save colony mask
            imageName = [wellFolder filesep 'colony' int2str(jj) '_colonyMask.tif'];
            imwrite(smoothMask, imageName);
        else
            finalImage = images{1};
        end
        figure; imshow(finalImage(:,:,dapiChannel),[]); title(['Colony' int2str(jj)]);
        figure; imshowpair(smoothMask, finalImage(:,:,dapiChannel)); title(['Colony' int2str(jj)]);
        
        
        % save colonyimages - all channels, dapi, nodal, lefty separately
        
        imageName = [wellFolder filesep 'colony' int2str(jj) '.tif'];
        % save multi channel images
        for nn = 1:size(finalImage,3)
            if nn == 1
                imwrite(finalImage(:,:,nn), imageName);
            else
                imwrite(finalImage(:,:,nn), imageName, 'WriteMode', 'append');
            end
            
        end
        
        % save dapi, spot channels
        for mm = [dapiChannel spotChannels membraneChannel]
            imageName = [wellFolder filesep 'colony' int2str(jj) '_ch' int2str(mm) '.tif'];
            imwrite(finalImage(:,:,mm), imageName);
        end
        
        
        colonyCounter = colonyCounter+1;
    end
    toc;
end

%%
%%
function finalImage = makeColonyImage(img, nImagesX, nImagesY, overlapLR, overlapTB)
% play with the overlaps to get the best match for your images
%% align images left-right
nImages = nImagesX*nImagesY;
rowImageIds = reshape([1:nImages], nImagesX, nImagesY)';
alignedImages = cell(1,size(rowImageIds,1));
side = 4;

alignmentChannel = 1;

for ii = 1:size(rowImageIds,1)

    img1 = img{rowImageIds(ii,1)};
    img2 = img{rowImageIds(ii,2)};
    alignedImage = alignTwoImagesFourier_multiChannel(img1, img2, side, overlapLR, alignmentChannel);
    
    if size(rowImageIds,2) == 3
        img3 = img{rowImageIds(ii,3)};
        alignedImages{ii} = alignTwoImagesFourier_multiChannel(alignedImage, img3,  side, overlapLR, alignmentChannel);
    else
        alignedImages{ii} = alignedImage;
    end
end
%% align images top-bottom
img1 = alignedImages{1};
img2 = alignedImages{2};

side = 1;
[img12, ~, ~] = alignTwoImagesFourier_multiChannel(img1, img2,  side, overlapTB, alignmentChannel);

if size(rowImageIds,1) == 3
    img1 = img12;
    img2 = alignedImages{3};
    [finalImage, ~, ~] = alignTwoImagesFourier_multiChannel(img1, img2,  side, overlapTB, alignmentChannel);
else
    finalImage = img12;
end
end
%%
