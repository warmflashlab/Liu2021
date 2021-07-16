close all
clear all
maxproj_dir='D:\2021-02-18-smFISH\Experiment2\MIP\E6_Wnt3a_2hr\';
save_dir = 'D:\2021-02-18-smFISH\Experiment2\';
test_str = '_E6Wnt3a_2h';%'Activin' num2str(v(tp)) 'h_MIP' 
chan_nm = {'Dapi','Lefty mRNA','Nodal mRNA'};
segm_chan = 1;% dapi
ff = dir(maxproj_dir);
q_chan =[2 3];% [1 2 3]
cd(maxproj_dir);
prob_thresh=0.9;
% params for kymograph for nodal movie
pxl_to_micron = 0.650;%( pxltomicron = 0.325 - 40X)
dilate = 5; % im pixels
small_stuff=200;

for ii=1:size(q_chan,2)
mask=[];
fnm = struct;
ilastik = struct;
mean_expr=[];
v = 1;
data = struct;
mask = struct;
q = 1;           
pos_dat=[];    
save_nm = [save_dir chan_nm{q_chan(ii)} '_E6Wnt3a_2h.mat'];

for jj=1:size(ff,1)
    if  ~isdir(ff(jj).name) &&  ~isempty(regexp(ff(jj).name,[test_str ],'ONCE')) && ~isempty(regexp(ff(jj).name,['_w000' num2str(segm_chan-1) '.tif'],'ONCE')) %'mTeSR'
        
        fnm(q).name =  ff(jj).name;
        ilastik(q).name = [fnm(q).name(1:end-4) '_Probabilities.h5'];% dapi
        %disp(ilastik(q).name);
        mask(q).all = imfill(bwareaopen(readIlastikProbMask(ilastik(q).name,prob_thresh),small_stuff),'holes');
        disp([ff(jj).name(1:size(ff(jj).name,2)-5) num2str(q_chan(ii)-1) '.tif']);
        img = imread([ff(jj).name(1:size(ff(jj).name,2)-5) num2str(q_chan(ii)-1) '.tif']);% stain
        img_no_bg = img-mean(img(~imdilate(mask(q).all,strel('disk',dilate))));
        figure,imshow(imdilate(mask(q).all,strel('disk',dilate)));
        %imshow(img_no_bg,[]);
        %imshowpair(img,imdilate(mask(q).all,strel('disk',dilate)));img,mask(q).all
        stats = regionprops(imdilate(mask(q).all,strel('disk',dilate)),img_no_bg,'Area','Centroid','MeanIntensity');
        size(stats)
        mean_expr = mean(cat(1,stats.MeanIntensity));
        data(q).expr = mean_expr;
        q = q+1;
    end
    
end
    chan = q_chan(ii);
    pos_dat = cat(1,data.expr);
    save(save_nm,'maxproj_dir','segm_chan','q_chan','v','pos_dat','chan_nm','small_stuff','prob_thresh','chan');

end
disp('done')

%% 
close all
ff = dir('D:\2021-02-18-smFISH\Experiment2');
test_str='Nodal';
w = 1;
for jj=1:size(ff,1)
    if  ~isdir(ff(jj).name) &&  ~isempty(regexp(ff(jj).name,[test_str ],'ONCE')) && ~isempty(regexp(ff(jj).name,'.mat','ONCE')) %'mTeSR'
        load(ff(jj).name);
        disp(ff(jj).name);
figure(1), errorbar(w,mean(pos_dat),std(pos_dat),'p','LineWidth',2);hold on

w = w+1;
    end
end
ylabel('Mean pixel intensity in dilated nuclei')
title(chan_nm{chan})
h6 = figure(1);
h6.CurrentAxes.FontSize = 12;
h6.CurrentAxes.LineWidth = 2;
xlim([0 10])
title(['smFISH ' test_str])    
legend('E6 CHIR 2h','E6 IWP2 2h','E6 Wnt3a 2h','E6 1hr','E6 2hr');

