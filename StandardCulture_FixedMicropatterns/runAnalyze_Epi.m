clear all; close all;

addpath(genpath('/Users/lizhongliu/Documents/MatlabStemcellsIdse')); 

dataDir = '/Volumes/Li/20200610Nodal/H336h';

% .vsi files
filelist = {'Process_567.vsi'};

vsifile = fullfile(dataDir,filelist{1});
[~,barefname,~] = fileparts(vsifile);

resultsDir = fullfile(dataDir,'results');
if ~exist(resultsDir,'dir')
    mkdir(resultsDir);
end
previewDir = fullfile(dataDir,'preview');
if ~exist(previewDir,'dir')
    mkdir(previewDir);
end
if ~exist(fullfile(resultsDir,'figs'),'dir')
    mkdir(fullfile(resultsDir,'figs'));
end

% TODO : local normalization by DAPI

% metadata
%----------------

meta = Metadata(vsifile);

% manually entered metadata
%------------------------------

meta.channelLabel = {'DAPI','NodalCit','pSmad2','NodalRNA'}
                     

meta.conditions = {'a'};

% IMPORTANT: we analyze one set at a time: change here
%antibodyset = 1;
%NA = numel(meta.channelLabel);
%filelist = filelist(antibodyset:NA:end);
%meta.channelLabel = meta.channelLabel{antibodyset};

% subset 
%conditionset = 1;
%filelist = filelist(conditionset:2:end);
%meta.conditions = meta.conditions(conditionset:2:end);

nucChannel = 1;
dataChannels = [1 2 3 4];

markerChannels = setdiff(dataChannels, nucChannel);

%% load data for parameter check

n = 1;
m = 2;
I = ones([1 numel(meta.conditions)]);
ymin = I*n*2048;
xmin = I*n*2048;
ymax = ymin + m*2048;
xmax = xmin + m*2048;
preview = {};
Ilim = {};

series = 1;

% load data and determine intensity limits per file
for fi = 1:numel(meta.conditions)

    vsifile = fullfile(dataDir,filelist{fi});
    metatmp = Metadata(vsifile);
    %channelOrder = orderChannels(metatmp)
    channelOrder = [1 2 3 4];
    
    xmax(fi) = min(xmax(fi), metatmp.xSize); 
    ymax(fi) = min(ymax(fi), metatmp.ySize);
    W = xmax(fi)-xmin(fi)+1; H = ymax(fi)-ymin(fi)+1;
    img_bf = bfopen_mod(vsifile,xmin(fi),ymin(fi),W,H,series);

    for ci = 1:numel(dataChannels)

        im = img_bf{1}{channelOrder(dataChannels(ci)),1};
        preview{ci, fi} = im;
        maxim = double(max(im(:)));
        minim = double(min(im(:)));
        tolerance = 0.04;
        Ilim{dataChannels(ci), fi} = stretchlim(mat2gray(im), tolerance)*(maxim-minim) + minim;
    end
end

% combined limits
for ci = 1:numel(dataChannels)  

    Imin = min(min([Ilim{dataChannels(ci),:}]));
    Imax = max(max([Ilim{dataChannels(ci),:}]));
    allIlim{dataChannels(ci)} = [Imin Imax]; %esto pone en min y max todo !!!!!
end

save(fullfile(resultsDir,['overviewIlim_' meta.channelLabel{dataChannels}]), 'Ilim', 'allIlim');

%% combine lookup tables and save RGB previews

%allIlim{1} = [400 500];
%allIlim{2} = [150 4000];
%allIlim{3} = [100 700];
%allIlim{4} = [100 500];

% use Ilim from other set
%bla = load(fullfile(resultsDir,['overviewIlim_' meta.channelLabel{dataChannels}]));
%allIlim = bla.Ilim;

preview8bit = {}; 
for fi = 1:numel(meta.conditions)

    for ci = 1:numel(dataChannels)  
        
        preview8bit{ci,fi} = uint8((2^8-1)*mat2gray(preview{ci, fi}, allIlim{dataChannels(ci)}));
    end

    previewChannels = 2:4;
    [~,barefname,~] = fileparts(filelist{fi});
    filename = fullfile(dataDir, 'preview', ['RGBpreview_' barefname '_' [meta.channelLabel{dataChannels}] '.tif']);
    RGBim = cat(3,preview8bit{previewChannels,fi});
    imwrite(RGBim, filename);
end

%% smaller previews to copy to slide

N = numel(meta.conditions);
n = ceil(N/4);
%m = ceil(N/n);
m = min(N,4);
screensize = get( 0, 'Screensize' );
margin = 50;
fs = 20;
w = screensize(3);
h = n*(screensize(3)/m + margin/2);
% if h > (screensize(4)-100)
%     w = w*(screensize(4)-100)/h;
%     h = screensize(4);
% end


figure('Position', [1, 1, w, h]);
for i = 1:2
    for fi = 1:N

        subplot_tight(n,m,fi)
        RGBim = cat(3,preview8bit{previewChannels,fi}); %maybe I need to change the order for different color
        %RGBnuc = repmat(preview8bit{dataChannels==nucChannel,fi},[1 1 3]); %preview8bit tiene dividida las imagenes 
        %RGBim = RGBim + 0.5*RGBnuc;
        %ims = {RGBim, RGBnuc};
        %imshow(ims{i});
        imshow(RGBim);
        titlestr = [meta.conditions{fi} ' (' filelist{fi}(end-7:end-4) ')'];
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none');
        labelstr = ['\color{red}'   meta.channelLabel{previewChannels(1)}...
                    '\color{green}' meta.channelLabel{previewChannels(2)}...
                    '\color[rgb]{.1, 0.5, 1}' meta.channelLabel{previewChannels(3)}]; %Magenta RGB 229, 9, 127 = [1 0 1]
        text(margin, size(RGBim,1) - 2.5*margin, labelstr,'FontSize',fs,'FontWeight','bold');
    end
    if i == 1
        fname = fullfile(resultsDir,['overview' meta.channelLabel{markerChannels} '.png']);
    else
        fname = fullfile(resultsDir,['overview' meta.channelLabel{markerChannels} 'DAPI.png']);
    end
    saveas(gcf, fname);
end
%close;
%%
% For single channel overviews:

singleChannels = dataChannels;

for currentChannel = singleChannels
    figure('Position', [1, 1, w, h]);
    for fi = 1:N

        subplot_tight(n,m,fi)
        RGBim = repmat(preview8bit{currentChannel,fi},[1 1 3]);
        imshow(RGBim);
        titlestr = meta.conditions{fi};
        title(titlestr,'FontSize',fs,'FontWeight','bold','Interpreter','none');
        labelstr = ['\color{red}'   meta.channelLabel{currentChannel}];
        text(margin, size(RGBim,1) - 2.5*margin, labelstr,'FontSize',fs,'FontWeight','bold');
    end
    fname = fullfile(resultsDir,['overview' meta.channelLabel{currentChannel} '.png']);
    saveas(gcf, fname);
end

%% load segmentation to test parameters %Ilastik needed

fi = 1;
vsifile = fullfile(dataDir, filelist{fi});

P = Position(meta.nChannels, vsifile, meta.nTime);
P.setID(pi);
dataSubdir = fileparts(vsifile);
seg = P.loadSegmentation(fullfile(dataSubdir,'MIP'), nucChannel);
segpart = seg(ymin(fi):ymax(fi), xmin(fi):xmax(fi));

%% finding the right options for extract data

% externally I will have all indices starting at 1
% the Andor offset to start at 0 will be internal

%blue is ilastik's and that has been discarded, pink is the one that stays.


opts = struct(  'cytoplasmicLevels',    false,... 
                    'dataChannels',     1:meta.nChannels,... 
                    'tMax',             meta.nTime,...
                    'segmentationDir',  fullfile(dataDir,'MIP'),...
                    'nucShrinkage',     2,...
                    'cytoSize',         3,...
                    'cytoMargin',       3,...
                    'bgMargin',         3);
                
opts.cleanupOptions = struct('separateFused', true,'openSize',4,...
    'clearBorder',true,'minAreaStd', 1, 'erodeSize',1,'minSolidity',0.9,...
    'minArea',5); % the value of min area has to be divisible by 5, or it won't work 

% try out the nuclear cleanup settings on some frame:
segpartclean = nuclearCleanup(segpart, opts.cleanupOptions);

% visualize
s = 0.2;
nucim = mat2gray(preview8bit{dataChannels==nucChannel,fi});

RGBsegCheck = cat(3,   nucim + s*mat2gray(segpartclean),...
                        nucim,...
                        nucim + s*segpart);
figure, imshow(RGBsegCheck)

filename = fullfile(resultsDir, [barefname 'cleansegCheck.tif']);
imwrite(RGBsegCheck, filename);

%% if not good, try dirty segmentation
%if previous seg is good, skip this


dirtyopts = struct();
% 3.5 microns is a good scale for the filters
dirtyopts.s = round(5/meta.xres);
dirtyopts.areacutoff = round(dirtyopts.s^2);
dirtyopts.mask = segpart; %dirtyopts.absolutethresh = 200;
dirtyopts.toobigNstd = 3;
nucseg = dirtyNuclearSegmentation(preview{dataChannels==nucChannel,fi}, dirtyopts);

% visualize nuclear mask
newmaskedge = nucseg - imerode(nucseg,strel('disk',1));
impp = imadjust(mat2gray(nucim));
overlay = cat(3, impp, impp+newmaskedge, impp+segpartclean);
imshow(overlay)

filename = fullfile(resultsDir, [barefname 'dirtysegCheck.tif']);
imwrite(overlay, filename);

opts.dirtyOptions = dirtyopts;

%% process the images

for fi = 1:numel(filelist)
    
    fi 
    vsifile = fullfile(dataDir, filelist{fi});
    [dataSubdir, barefname,~] = fileparts(vsifile);
    metatmp = Metadata(vsifile);
    channelOrder = [1 2 3 4];%orderChannels(metatmp);
    opts.dataChannels = intersect(channelOrder, dataChannels,'stable');

    P = Position(meta.nChannels, vsifile, meta.nTime);
    P.setID(fi);
    
    opts.segmentationDir = fullfile(dataSubdir,'MIP');
    
    debugInfo = P.extractData(dataSubdir, channelOrder(nucChannel), opts);

    % save segcheck im
    %im = mat2gray(preview8bit{dataChannels==nucChannel,fi});
    im = cat(3,preview8bit{previewChannels,fi});
    mask = debugInfo.nucmask(ymin(fi):ymax(fi),xmin(fi):xmax(fi));
    mask = imdilate(mask,strel('disk',3)) - mask > 0;
    im(repmat(mask,[1 1 3])) = 255;
    filename = fullfile(resultsDir, [barefname 'segCheck.tif']);
    imwrite(im, filename);

    save(fullfile(resultsDir, [barefname '_positions.mat']),'P');
end

%% scatter cell locations on top of part of preview to check result

XY = P.cellData.XY;

inpreview =     XY(:,1) < xmax(fi) & XY(:,1) > xmin(fi)...
                & XY(:,2) < ymax(fi) & XY(:,2) > ymin(fi);
XY = XY(inpreview,:);

nucim = mat2gray(preview8bit{dataChannels==nucChannel,fi});
imshow(nucim,[])
hold on
scatter(XY(:,1) - xmin(fi), XY(:,2) - ymin(fi), 200, '.r');
hold off


%% make distributions

distChannels = 1:4;

% load all data  
allData = {};
for fi = 1:numel(filelist)

    vsifile = fullfile(dataDir, filelist{fi});
    [~,barefname,~] = fileparts(vsifile);
    load(fullfile(dataDir,'results', [barefname '_positions.mat']));
    allData{fi} = P;
end

% make distributions out of processed data
stats = cellStats(allData, meta, distChannels);
tolerance = 0.01;
nbins = 50;
stats.makeHistograms(nbins, tolerance);

%% plot distributions

N = numel(distChannels);
n = ceil(N/4); m = min(N,4);
margin = 10;
screensize = get( 0, 'Screensize' );
h1 = figure;
h2 = figure('Position', [1, 1, screensize(3), n*(screensize(3)/m + margin/2)]);

% overlay of the distribution of different conditions
for channelIndex = distChannels

    figure(h1)
    clf
    stats.plotDistributionComparison(channelIndex)
    ylim([0 0.35]);

    filename = ['distOverlay_' meta.channelLabel{channelIndex}];
    saveas(gcf, fullfile(resultsDir, 'figs', filename));
    saveas(gcf, fullfile(resultsDir, [filename '.png']));

    clf
    cytoplasmic = false;
    cumulative = true;
    stats.plotDistributionComparison(channelIndex, cytoplasmic, cumulative)
    
    filename = ['cumdistOverlay_' meta.channelLabel{channelIndex}];
    saveas(gcf, fullfile(resultsDir, 'figs', filename));
    saveas(gcf, fullfile(resultsDir, [filename '.png']));
    
    figure(h2)
    subplot_tight(n,m, channelIndex, [3 1]*0.05)
    stats.plotDistributionComparison(channelIndex)
    ylim([0 0.35]);
end

filename = ['distoverlayAll_' stats.channelLabel{distChannels}];
saveas(figure(h2),fullfile(resultsDir,'figs', filename));
saveas(figure(h2),fullfile(resultsDir, [filename '.png']));

%% try clustering in 2D

k = 2; %how many clusters you wanna find
regularizationValue = 10^(-3);
cutoff = 0.8;
[y, x] = stats.makeClusters(k, regularizationValue, cutoff);
% ci1 = 1; ci2 = 2;
% scatter(x(:,ci1),x(:,ci2),1,y)

%% 

for ci1 = 1:numel(distChannels)
    for ci2 = ci1+1:numel(distChannels)

        figure,
        c1 = distChannels(ci1);
        c2 = distChannels(ci2);
        showClusters = true;
        for conditionIdx = 1:numel(stats.conditions)
            stats.makeScatterPlot(conditionIdx, [c1 c2], showClusters)
        end
        title('all data clustered');

        [~,barefname,~] = fileparts(filelist{conditionIdx});
        filename = ['cluster_' meta.channelLabel{c1} '_' meta.channelLabel{c2}];
        saveas(gcf, fullfile(resultsDir, 'figs', filename));
        saveas(gcf, fullfile(resultsDir, [filename '.png']));
    end
end

%% if bad manual cluster definition
%you can skip this if you like with the previous one

k = 1;
ci1 = 2; ci2 = 3;
stats.makeClustersManual(k, ci1, ci2)

%%

load(fullfile(resultsDir,'stats'),'stats');

%% if good visualize clusters for each condition

N = numel(meta.conditions);
n = ceil(N/4); m = min(N,4);
margin = 10;
screensize = get( 0, 'Screensize' );
h1 = figure;
h2 = figure('Position', [1, 1, screensize(3), n*(screensize(3)/m + margin/2)]);

for ci1 = 1:numel(distChannels)
    for ci2 = ci1+1:numel(distChannels)
        
        channelIdx = distChannels([ci1 ci2]);
        showClusters = true; % SET TO TRUE IF CLUSTERS
        
        figure(h2)
        clf
        
        for conditionIdx = 1:numel(stats.conditions)

            figure(h1)
            clf
            stats.makeScatterPlot(conditionIdx, channelIdx, showClusters)
            
            [~,barefname,~] = fileparts(filelist{conditionIdx});
            filename = ['cluster_' stats.channelLabel{channelIdx(1)} '_'...
                        stats.channelLabel{channelIdx(2)} '_' barefname];
            saveas(gcf, fullfile(resultsDir, 'figs', filename));
            saveas(gcf, fullfile(resultsDir, [filename '.png']));
            
            figure(h2)
            subplot_tight(n,m, conditionIdx, [2 1]*0.05)
            stats.makeScatterPlot(conditionIdx, channelIdx, showClusters)
        end
        
        filename = ['clusterAll_' stats.channelLabel{channelIdx(1)} '_'...
                        stats.channelLabel{channelIdx(2)}];
        saveas(figure(h2),fullfile(resultsDir,'figs', filename));
        saveas(figure(h2),fullfile(resultsDir, [filename '.png']));
    end
end

%% this doesn't work
% plot distributions per cluster
stats.makeClusterHistograms();

for clusterLabel = 1:stats.nClusters
    for channelIndex = distChannels

        cumulative = false;
        stats.plotDistributionComparison(channelIndex, cumulative, clusterLabel)

        filename = ['distOverlay_' meta.channelLabel{channelIndex} '_cluster' num2str(clusterLabel)];
        saveas(gcf, fullfile(resultsDir, 'figs', filename));
        saveas(gcf, fullfile(resultsDir, [filename '.png']));
        close;

        cumulative = true;
        stats.plotDistributionComparison(channelIndex, cumulative, clusterLabel)

        filename = ['cumdistOverlay_' meta.channelLabel{channelIndex} '_cluster' num2str(clusterLabel)];
        saveas(gcf, fullfile(resultsDir, 'figs', filename));
        saveas(gcf, fullfile(resultsDir, [filename '.png']));
        close;
    end
end
%%
% summary 
txtfile = fullfile(resultsDir, ['statsSummary_' meta.channelLabel{dataChannels} '.txt']);
if exist(txtfile,'file')
    delete(txtfile);
end
diary(txtfile);
diary on
stats.printSummary();
diary off

% save stats object
save(fullfile(resultsDir,['stats_' meta.channelLabel{dataChannels}]),'stats');
