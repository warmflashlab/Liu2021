function [radialProfile_avg, colonyCounter] = computeCombinedRadialAverage(colonies, goodColoniesId, meta)
%% calculate avg radial profile using good colonies

%% ------------ Inputs ----------
% colonies: object of type colonyS
% goodColoniesId: position of good colonies in the colonies object

%% ----------- Output -----------
% radialProfile_avg: structure with combined radial average, standard
% deviation and standard error for both dapiNormalized and not normalized
% values.
% colonyCounter: number of colonies that contributed to each bin.

%%
if ~exist('goodColoniesId', 'var')
    goodColoniesId = 1:numel(colonies);
end

%% make rA_mean, rA_std, nPixels table
nBins = numel(colonies(2).radialProfile.bins)-1;
nChannels = size(colonies(1).radialProfile.dapiNormalized.mean,1);
membraneChannel = find(ismember(upper(meta.channelLabel), 'NONMEMBRANE'), 1);
rP = cell(2,2); % radialProfile [dapiNormalized_mean dapiNormalized_std; notNormalized_mean notNormalized_std]

nPixels = zeros(numel(goodColoniesId), nBins); % pixels table
mPixels = nPixels;
counter = 1; % colonyCounter

for ii = goodColoniesId
    for jj = 1:nChannels
        
        for dN = [1 2] %[dapiNormalized; notNormalized]
            if dN == 1
                rA = colonies(ii).radialProfile.dapiNormalized.mean(jj,:);
                std1 = colonies(ii).radialProfile.dapiNormalized.std(jj,:);
            else
                rA = colonies(ii).radialProfile.notNormalized.mean(jj,:);
                std1 = colonies(ii).radialProfile.notNormalized.std(jj,:);
            end
            
            if numel(rA)<nBins
                rA(nBins) = 0;  std1(nBins) = 0;  %filler - fills in zeros if bins are empty.
            end
            
            rP{dN,1}(counter,:,jj) = rA;
            rP{dN,2}(counter,:,jj) = std1;
            
        end
    end
    nPixels(counter,:) = colonies(ii).radialProfile.pixels;
    if ~isempty(membraneChannel)
        mPixels(counter,:) = colonies(ii).radialProfile.pixels;
    end
    counter = counter+1;
end

%%
colonyCounter = sum(nPixels~=0,1); % number of colonies

%% compute combined mean, std deviation
for ii = 1:2
    mean1 = sum(rP{ii,1}.*nPixels,1)./sum(nPixels,1);
    
    term1 = sum(rP{ii,2}.^2.*(nPixels-1) + rP{ii,1}.^2.*(nPixels),1);
    term2 = sum(nPixels,1).*(mean1.^2);
    std1 = sqrt((term1-term2)./(sum(nPixels,1) - 1));
    stdError1 = std1./sqrt(colonyCounter);
    
    mean1 = squeeze(mean1); mean1(any(isnan(mean1),2),:) = [];
    std1 = squeeze(std1); std1(any(isnan(std1),2),:) = [];
    stdError1 = squeeze(stdError1); stdError1(any(isnan(stdError1),2),:) = [];
    
    if ii == 1
        %dapi normalized
        radialProfile_avg.dapiNormalized.mean = mean1';
        radialProfile_avg.dapiNormalized.std = std1';
        radialProfile_avg.dapiNormalized.stdError = stdError1';
    else
        % not normalized
        radialProfile_avg.notNormalized.mean = mean1';
        radialProfile_avg.notNormalized.std = std1';
        radialProfile_avg.notNormalized.stdError = stdError1';
        
    end
end


if ~isempty(membraneChannel)
    mc = membraneChannel;
    mean1 = sum(rP{ii,1}(:,:,mc).*mPixels,1)./sum(mPixels,1);
    term1 = sum(rP{ii,2}(:,:,mc).^2.*(mPixels-1) + rP{ii,1}(:,:,mc).^2.*(mPixels),1);
    term2 = sum(mPixels,1).*(mean1.^2);
    std1 = sqrt((term1-term2)./(sum(mPixels,1) - 1));
    stdError1 = std1./sqrt(colonyCounter);
    
    mean1 = squeeze(mean1); mean1(any(isnan(mean1),2),:) = [];
    std1 = squeeze(std1); std1(any(isnan(std1),2),:) = [];
    stdError1 = squeeze(stdError1); stdError1(any(isnan(stdError1),2),:) = [];
    
    radialProfile_avg.notNormalized.mean(mc,:) = mean1';
    radialProfile_avg.notNormalized.std(mc,:) = std1';
    radialProfile_avg.notNormalized.stdError(mc,:) = stdError1';
    
    
end
end