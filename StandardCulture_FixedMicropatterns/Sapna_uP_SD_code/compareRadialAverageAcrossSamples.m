
%% ------------------------------- compare rA values across samples

%%
clearvars;
masterFolder = '/Volumes/Li/20200609/allData';
metadata  = [masterFolder filesep 'metadata.mat'];
load(metadata);
%%
samplesFolder = [masterFolder filesep 'processedData'];
samples = dir(samplesFolder);
toKeep = find(~cell2mat(cellfun(@(c)strcmp(c(1),'.'),{samples.name},'UniformOutput',false))); % remove non-named folders
samples = samples(toKeep);
[~, idx] = natsortfiles({samples.name});
samples = samples(idx);
%%
%load('/Users/sapnachhabra/Desktop/CellTrackercd/CellTracker60X/colorblind_colormap/goodColors_new.mat');
% 11 distinct colors. For more colors, use colorcube
colors = colorcube;
%%
% save rA values
%samples = samples([1 5 24]);
rA = cell(1,numel(samples));
nColonies = zeros(1,numel(samples));

for ii = 1:numel(samples)
    ii;
    %for ii = 1:numel(samples)
    %for ii = 1:numel(samples)
    outputFile = [samplesFolder filesep samples(ii).name filesep 'output.mat'];
    load(outputFile, 'radialProfile_avg', 'goodColoniesId', 'xValues');
    nColonies(ii) = numel(goodColoniesId);
    rA{ii} = radialProfile_avg;
end

%%
% plot
controlSample = 8;
dapiEqualize = 0; % make dapi max the same across conditions

%colors = cell2mat({[0.7 0 0.5]; [0 0.7 0];  [0 0 0.7]; [0.7 0 0.5]; [0 0.7 0];  [0 0 0.7]});
%colors = cell2mat({[1.0 0.6 0]; [1 0 0.9];  [0.7 0.7 0.7];  [0.7 0.5 0.0]; [0 0.5 0.5]});
%sampleLabels = {'dish:media 45h' 'MP:bmp 45h', 'MP:media 11h', 'MP:media 45h'};

sampleLabels = strrep({samples.name}, '_', ':');

if dapiEqualize == 1
    rA1_dapi = rA{controlSample}.notNormalized.mean(1,:)./max(rA{controlSample}.notNormalized.mean(1,:));
    rA_max = max(rA{controlSample}.notNormalized.mean./rA1_dapi,[],2);
else
    rA_max = max(rA{controlSample}.dapiNormalized.mean,[],2);
end
%%
sampleLabels = strrep({samples.name}, '_', ':');
%sampleLabels = {'bmp4', 'Dkk1:300ng/ml', 'Dkk1:2:300ng/ml', 'Dkk1:600ng/ml', 'Dkk1:900ng/ml', 'Lefty:500ng/ml', 'Lefty:750ng/ml',  'Lefty:1000ng/ml'};
%%
%% ------------- one channel - all samples
samplesToPlot = [1 2 3 4 5 6 7]; channelsToPlot = [1:4];
% plot
for ii = channelsToPlot
    figure; hold on; title(meta.channelNames{ii});
    
    %make legends of the right color.
    for jj = 1:numel(samplesToPlot)
        plot(0, 0, 'Color', colors(samplesToPlot(jj),:)); % font color
    end
    
    [~,hObj]=legend(sampleLabels(samplesToPlot), 'FontSize', 25, 'FontWeight', 'bold');   % return the handles array
    hL=findobj(hObj,'type','line');  % get the lines, not text
    set(hL,'linewidth',6);
    hT=findobj(hObj,'type','text');
    
    for jj = 1:numel(samplesToPlot)
        set(hT(jj),'Color', colors(samplesToPlot(jj),:)); % line color
    end
    legend('boxoff');
    
    counter = 1;
    for jj = samplesToPlot
        if dapiEqualize == 1
            rA1_dapi = rA{jj}.notNormalized.mean(1,:)./max(rA{jj}.notNormalized.mean(1,:));
            rA1 = rA{jj}.notNormalized.mean(ii,:)./rA1_dapi;
            rA1 = rA1./rA_max(ii);
            stdError1 = rA{jj}.notNormalized.stdError(ii,:)./rA1_dapi;
            stdError1 = stdError1./rA_max(ii);
        else
            rA1 = rA{jj}.dapiNormalized.mean(ii,:)./rA_max(ii);
            stdError1 = rA{jj}.dapiNormalized.stdError(ii,:)./rA_max(ii);
        end
        
        nBins = size(rA1,2);
        plot(xValues(1:nBins), rA1, 'Color', colors(jj,:), 'LineWidth', 5);
        errorbar(xValues(1:nBins), rA1, stdError1, 'Color', colors(jj,:), 'LineWidth', 2);
        xlabel('Distance from edge (\mum)'); ylabel('Intensity (a.u.)');
        ylim([0.2 3]);
        ax = gca; ax.FontSize = 25; ax.FontWeight = 'bold';
    end
end
%%
%% ----------------------------- one sample - all channels --------------------------------------
%%
dapiChannel = find(cell2mat(cellfun(@(c)strcmp(c,'DAPI'),upper(meta.channelNames),'UniformOutput',false)));
%% ----------------------------------------------------------------------------------------------
%%
samplesToPlot = [8:14]; channelsToPlot = [2:4];
for ii = samplesToPlot
    figure; hold on; title(samples(ii).name);
    
    for jj = 1:numel(channelsToPlot)
        plot(0,0,'Color', colors(channelsToPlot(jj),:));
    end
    
    [~,hObj]=legend(meta.channelNames(channelsToPlot), 'FontSize', 25, 'FontWeight', 'bold');   % return the handles array
    hL=findobj(hObj,'type','line');  % get the lines, not text
    set(hL,'linewidth',6);
    hT=findobj(hObj,'type','text');
    
    for jj = 1:numel(channelsToPlot)
        set(hT(jj),'Color', colors(channelsToPlot(jj),:)); % line color
    end
    legend('boxoff');
    
    counter = 1;
    for jj = channelsToPlot
        
        rA_max = max(rA{ii}.dapiNormalized.mean(jj,:));
        rA1 = rA{ii}.dapiNormalized.mean(jj,:)./rA_max;
        stdError1 = rA{ii}.dapiNormalized.stdError(jj,:)./rA_max;
        
        nBins = size(rA1,2);
        plot(xValues(1:nBins), rA1, 'Color', colors(jj,:), 'LineWidth', 5);
        errorbar(xValues(1:nBins), rA1, stdError1, 'Color', colors(jj,:), 'LineWidth', 2);
        xlabel('Distance from edge (\mum)'); ylabel('Intensity (a.u.)');
        ylim([0 1.5]);
        ax = gca; ax.FontSize = 25; ax.FontWeight = 'bold';
    end
    
end
%%
saveInPath =[masterFolder filesep 'rA_plots'];
saveAllOpenFigures(saveInPath);
%% ---------------------------------------------------------------------------------------------
%% ---------------------------- half max - front, back positions -------------------------------
%% ---------------------------------------------------------------------------------------------
samplesToPlot = 1:3;
counter1 = 1;
rA1_halfMax_back = zeros(numel(channelsToPlot), numel(samplesToPlot));
rA1_halfMax_front = rA1_halfMax_back;

for ii = channelsToPlot
    counter2 = 1;
    for jj = samplesToPlot
        rA1 = rA{jj}.dapiNormalized.mean(ii,:);
        rA1_halfMax = min(rA1) + (max(rA1) - min(rA1))/2;
        L1 = [xValues(1:nBins); rA1(1:nBins)];
        L2 = [xValues(1:nBins); repmat(rA1_halfMax, 1, nBins)];
        P = InterX(L1, L2);
        rA1_halfMax_back(counter1, counter2) = P(1,1);
        if size(P,2) == 1
            rA1_halfMax_front(counter1, counter2) = P(1,1);
        else
            rA1_halfMax_front(counter1, counter2) = P(1,2);
        end
        counter2 = counter2+1;
    end
    counter1 = counter1+1;
end

%%
rA1_halfMax_back(3,3) = 5;

%%
% Nodal, Lefty, smad2 front end movements

%%
figure; 
hold on;
x1 = [24 30 36];
for ii = 1:3
    y1 = rA1_halfMax_front(ii,:);
    plot(x1, y1, '*-', 'LineWidth', 5, 'Color', colors(ii+1,:));
end
xlabel('Time(h)'); ylabel('HalfMax:Position(\mum)');
legend(meta.channelNames(channelsToPlot));
ax = gca;
ax.FontSize = 30;
ax.FontWeight = 'bold';
%%
for ii = 1:3
    y1 = rA1_halfMax_back(ii,:);
    plot(x1, y1, '*-', 'LineWidth', 3, 'Color', colors(ii+1,:));
end
%%
% width of half_max region
figure; hold on;
for ii = 1:3
    y1 = rA1_halfMax_front(ii,:) - rA1_halfMax_back(ii,:);
    plot(x1, y1, '*-', 'LineWidth', 5, 'Color', colors(ii+1,:));
end
xlabel('Time(h)'); ylabel('HalfMax:Width (\mum)');
legend(meta.channelNames(channelsToPlot));
ax = gca;
ax.FontSize = 30;
ax.FontWeight = 'bold';

%%
saveInPath =[masterFolder filesep 'halfMaxPositions'];
saveAllOpenFigures(saveInPath);
%%







