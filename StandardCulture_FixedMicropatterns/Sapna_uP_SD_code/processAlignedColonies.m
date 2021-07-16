
%% ---------------------------------------------------------------------------
% saves radialaverage information for all colonies
close all;
clearvars;
addpath(genpath('/Users/lizhongliu/Documents/MatlabFunctions'));
masterFolder =  '/Volumes/Li/20200609/allData';
tooHighIntensity = 15000; % cut-off for removing super bright non-nuclear pixels, LSM (4096), epi (8500)

metadata  = [masterFolder filesep 'metadata.mat'];
load(metadata);

samplesFolder = [masterFolder filesep 'processedData'];
samples = dir(samplesFolder);

toKeep = find(~cell2mat(cellfun(@(c)strcmp(c(1),'.'),{samples.name},'UniformOutput',false))); % remove non-named folders
samples = samples(toKeep);
[~, idx] = natsortfiles({samples.name});
samples = samples(idx);
%%
dapiChannel = find(cell2mat(cellfun(@(c)strcmp(c,'DAPI'),upper(meta.channelNames),'UniformOutput',false)));

outerBin = 10; % in microns
bins = getBinEdgesConstantArea(meta.colRadiiMicron(1), outerBin); % bin area is constant.

% for live cell imaging, specify imaging break points
%cMask_timepoints = [1 12];
%%
% --- save output files for each sample
%for ii = 1
%for ii = 5
for ii = 1:size(samples,1)
    tic;
    disp(['Sample ' int2str(ii)]);
    
    files = dir([samples(ii).folder filesep samples(ii).name filesep '*_ch' int2str(dapiChannel) '_Simple Segmentation.h5']); % gives unique samples
    nFiles = size(files,1);
    colonies = colonyS;
    
    %for jj = 5
    for jj = 1:nFiles
        colonies(jj) = colonyS(jj,files(1).folder, ['colony' int2str(jj) '.tif'], meta.colRadiiMicron);
        colonies(jj) = calculateRadialProfile(colonies(jj), meta, bins, tooHighIntensity); %
        %fixed cell
        %colonies(jj) = calculateRadialProfileD(colonies(jj), meta, bins, cMask_timepoints); %live cell
    end
    
    outputFile = [files(1).folder filesep 'output.mat'];
    save(outputFile, 'colonies');
    toc;
end

%%
%% ---------------------------------------------------------------------------------
%% ---------------------------------------------------------------------------------
%% ----------------- filter ---------------------
% filter colonies, plot radial average profile
%for sampleId = 1:size(sa
%for sampleId = 9:numel(samples)
sampleId = 1;
clearvars -except meta samples samplesFolder sampleId
close all;

outputFile = [samplesFolder filesep samples(sampleId).name filesep 'output.mat'];
load(outputFile);

%% --------------------------------------------------------------------------------
%% --------------------------------------------------------------------------------
%% ---------------- plot -----------------------
% plot
% 1) indivdual colonies - ch1,2,3
% 2) average

%% -- plot individual colony intensity profiles
% (fixed cell data)

% 1)
nColonies = numel(colonies);
bins = colonies(1).radialProfile.bins;
xValues = (bins(1:end-1)+bins(2:end))/2;

figure;
% non-normalized values
for ii = 1:meta.nChannels
    subplot(2,2,ii); hold on;
    for jj = 1:nColonies
        rA = colonies(jj).radialProfile.notNormalized.mean(ii,:);
        if jj < 6
            lw = 1;
        elseif jj < 12
            lw = 4;
        else
            lw = 7;
        end
        plot(xValues(1:numel(rA)), rA, 'LineWidth', lw);
    end
    legend(strcat('Colony',strsplit(int2str(1:nColonies), ' '))); title(meta.channelNames{ii});
end

figure;
% normalized values
coloniesToPlot = 1:nColonies;
for ii = 1:meta.nChannels
    subplot(2,2,ii); hold on;
    for jj = coloniesToPlot
        rA = colonies(jj).radialProfile.dapiNormalized.mean(ii,:);
        if jj <6
            lw = 1;
        elseif jj < 12
            lw = 4;
        else
            lw = 7;
        end
        plot(xValues(1:numel(rA)), rA, 'LineWidth', lw);
    end
    legend(strcat('Colony',strsplit(int2str(coloniesToPlot), ' '))); title(meta.channelNames{ii});
end
%%
% filter colonies
badColonies = []; % for no badcolonies, set badColonies = [];
%%
% 2) calculate combined radial profile using good colonies
% make mean, std, nPixels table
if ~isempty(badColonies)
    goodColoniesId = setxor([1:numel(colonies)], badColonies);
else
    goodColoniesId = 1:numel(colonies);
end


radialProfile_avg = computeCombinedRadialAverage(colonies, goodColoniesId, meta);
%
%
%  ---- plot radialavg
% std deviation
% dapiNormalized
xValues = 0.5*(colonies(1).radialProfile.bins(1:end-1)+colonies(1).radialProfile.bins(2:end));
nValues = size(radialProfile_avg.dapiNormalized.mean,2);
figure;
for ii = 1:meta.nChannels
    subplot(2,2,ii);
    errorbar((1:nValues), radialProfile_avg.dapiNormalized.mean(ii,:), radialProfile_avg.dapiNormalized.std(ii,:));
    title(meta.channelNames{ii});
end

% notNormalized
xValues = 0.5*(colonies(1).radialProfile.bins(1:end-1)+colonies(1).radialProfile.bins(2:end));
figure;
for ii = 1:meta.nChannels
    subplot(2,2,ii);
    errorbar(xValues(1:nValues), radialProfile_avg.notNormalized.mean(ii,:), radialProfile_avg.notNormalized.std(ii,:));
    title(meta.channelNames{ii});
end


% std error
% dapiNormalized
figure;
for ii = 1:meta.nChannels
    subplot(2,2,ii);
    errorbar(xValues(1:nValues), radialProfile_avg.dapiNormalized.mean(ii,:), radialProfile_avg.dapiNormalized.stdError(ii,:));
    title(meta.channelNames{ii});
end

% notNormalized
xValues = 0.5*(colonies(1).radialProfile.bins(1:end-1)+colonies(1).radialProfile.bins(2:end));
figure;
for ii = 1:meta.nChannels
    subplot(2,2,ii);
    errorbar(xValues(1:nValues), radialProfile_avg.notNormalized.mean(ii,:), radialProfile_avg.notNormalized.stdError(ii,:));
    title(meta.channelNames{ii});
end

% std error dapi normalized, all channels
figure; hold on;
channelIds = 2:meta.nChannels;
for ii = channelIds
    rA = radialProfile_avg.dapiNormalized.mean(ii,:);
    rA_max = max(rA);
    rA1 = rA/rA_max;
    rAerror1 = radialProfile_avg.dapiNormalized.stdError(ii,:)./rA_max;
    errorbar(xValues(1:nValues), rA1, rAerror1);
end
xlabel('Distance from edge (mum)'); ylabel('Intensity (a.u.)');
legend(meta.channelNames(channelIds));

% save
save(outputFile, 'radialProfile_avg', 'goodColoniesId', 'xValues', '-append');
saveInPath = [samplesFolder filesep samples(sampleId).name filesep 'rA_plots'];
saveAllOpenFigures(saveInPath);

%end
%%



























