%example run
addpath(genpath('/Users/ll/Documents/MATLAB/MatlabFunctions')); 
p = runMovieSimpleSegmentation('MaxMemExp_G001.tif',20,100,1:2);
%%
outNuc = [p.pixAvgData.nucAvg];
outNucCyt = [p.pixAvgData.nucCytAvg];

figure; 
%plot(outNuc(1,:),'r.-'); hold on;
plot(outNucCyt(2,:)./outNuc(1,:),'b.-'); hold off;
xlabel('time (frames)'); ylabel('Intensity (a.u)');
%%
tic;
basefile = 'MaxMemExp';
for ii = 1:32
    disp(int2str(ii));
    filename = sprintf('%s_G%.3d.tif',basefile,ii);
    p = runMovieSimpleSegmentation(filename,20,100,1:2);
    positions(ii) = p;
end
toc;
%% non-normalized to control condition (mTeSR)
% filenrs = {1:4 ,5:8,7:12,13:16,17:20,21:24,25:28,29:32,33:36};
% conds = {'w0-a0','w0-a0.075','w0-a0.15','w0-a0.3', 'w0-a0.6','w0-a1.25','w0-a2.5','w0-a5','w0-a10'};
% figure; hold on;
% for ii = 1:length(filenrs)
%     pos_now = positions(filenrs{ii});
%     for jj = 1:length(pos_now)
%         tt = [pos_now(jj).pixAvgData.nucCytAvg];
%         tt2 = [pos_now(jj).pixAvgData.nucAvg];
%         allavgs(jj,:) = tt(2,:);
%     end
% 
%     avg = mean(allavgs);
%     err = std(allavgs);
%     
%     rr = 1:length(avg)*35;
%     errorbar(rr,avg,err,'LineWidth',2);
% end
% xlabel('time (min)');
% ylabel('intensity (au)');
% legend(conds,'Location','Best');
% set(gca,'FontSize',18);
% xlim([1 27]);
% saveas(gcf,'NodalQuant.png');

%% Normalized to control condition
filenrs = {1:4 ,5:8,7:12,13:16,17:20,21:24,25:28,29:32};
conds = {'mTeSR','A10','CHIR24h','CHIR5h','Wnt24','Wnt5h','A10+CHIR24','A10+Wnt'};
figure; hold on;
controlcondition = 1;
for ii = 1:length(filenrs)
    pos_now = positions(filenrs{ii});
    for jj = 1:length(pos_now)
        tt = [pos_now(jj).pixAvgData.nucCytAvg];
        tt2 = [pos_now(jj).pixAvgData.nucAvg];
        allavgs(jj,:) = tt(2,:);
    end

    avg = mean(allavgs);
    err = std(allavgs)/sqrt(length(allavgs)); %SEM
%   err = std(allavgs); %SD

    if ii==controlcondition
        avgcontrol = avg;
    end

    timeinterval = 35; %minutes
    rr = (1:length(avg))*timeinterval/60;
   
     errorbar(rr,avg./avgcontrol,err./(avgcontrol.^1),'LineWidth',2);
     


end

xlabel('Time (h)');
ylabel('Normalized intensity');
legend(conds,'Location','Best');
set(gca,'FontSize',24,'Box',1,'LineWidth',3);
xlim([1 30]);
ylim([0.7 2.5]);
saveas(gcf,'NodalQuant.png');
save('data.mat');
