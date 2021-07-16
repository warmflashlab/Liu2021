%% get the number of mRNA spots as a function of distance from the border between receiving and producing cells

%1. get the meberane mask 
close all
clear all
maxproj_dir= 'D:\2021-01-01-smFISH\Dapi_Lefty_Produce_Nodal\MIP\WT_jux_WT_wnti\20hr\';
spots_dir='D:\2021-01-01-smFISH\Dapi_Lefty_Produce_Nodal\MIP\WT_jux_WT_wnti\Spots_20hr\'; 
save_nm = 'D:\2021-01-01-smFISH\Dapi_Lefty_Produce_Nodal\MIP\WT_jux_WT_wnti\WT_jux_WT_wnti_20hr.mat';
% load the .csv file manually (for now) and save into mat file
%   pos = 4;
%   save([spots_dir 'p000' num2str(pos-1) '_spots.mat'],'p0003w0003');

ff = dir(maxproj_dir);
segm_chan = 2;        % channel used to get the segementation (masks of receiving/producing cells areas)
receivers_label = 0;  % which cells are labelled with membrane marker
chan_nm = 'Nodal_mRNA';
cd(maxproj_dir);
mask_producing = [];
mask_receiving= [];
fnm = struct;
ilastik = struct;
prob_thresh=0.1;
% params for kymograph for nodal movie
pxl_to_micron = 0.650;%( pxltomicron = 0.325 - 40X)
dxy = 5;% size of slices to look at incrementally, away from the borderin microns
max_d = 200;%  max distance from the border ,in microns
small_stuff=200;
mRNA_perArea=struct;
v = [1 2 3 4];
img = struct;
q = 1;

for pos =1:size(v,2)%1:size(positions,2)    
for jj=1:size(ff,1)
    if ~isdir(ff(jj).name) &&  ~isempty(regexp(ff(jj).name,['p000' num2str(pos-1) '_w000' num2str(segm_chan) ],'ONCE')) && ~isempty(strfind(ff(jj).name,'.tif'))
        fnm(q).name = ff(jj).name;% this is the image/movie name
        ilastik(q).name = [fnm(q).name(1:end-4) '_Probabilities.h5'];
        img(q).im = imread(fnm(q).name);
        disp(ilastik(q).name);
        if receivers_label
        mask_receiving = imfill(bwareaopen(readIlastikProbMask(ilastik(q).name,prob_thresh),small_stuff),'holes');
        mask_producing =  bwareaopen(~mask_receiving,small_stuff); 
        else
        mask_producing = imfill(bwareaopen(readIlastikProbMask(ilastik(q).name,prob_thresh),small_stuff),'holes');
        mask_receiving =  bwareaopen(~mask_producing,small_stuff);     
            
        end
        q = q+1;
        
    end      
end
stats = regionprops(mask_producing,'Area');
mask_producing =  imdilate(bwareaopen(mask_producing,(min(cat(1,stats.Area))+10)),strel('disk',5)); 
mask_receiving =  ~mask_producing;  
figure,imshowpair(mask_receiving,mask_producing)


if ~isempty(intersect(pos-1,v-1))
spots_data=importdata([spots_dir 'p000' num2str(pos-1) '_spots.mat']);
['p000' num2str(pos-1) '_spots.mat']
spots_data(spots_data==0)=1;
img_mrna = zeros(size(mask_receiving));
for jj=1:size(spots_data,1)
img_mrna(spots_data(jj,1),spots_data(jj,2))=1;
end
end
imshow(img_mrna)
hold on
plot(spots_data(:,2),spots_data(:,1),'p','LineWidth',0.05);hold on

[expr_dat_inreceiving,dist_i]=get_mrna_dist(pxl_to_micron,mask_producing,dxy,max_d,img_mrna);
mRNA_perArea(pos).rcv = [dist_i(1:size(expr_dat_inreceiving,1))'   expr_dat_inreceiving ];

[expr_dat_inprod,dist_i]=get_mrna_dist(pxl_to_micron,mask_receiving,dxy,max_d,img_mrna);
mRNA_perArea(pos).prdc = [dist_i(1:size(expr_dat_inprod,1))'   expr_dat_inprod ];
%close all
end

save(save_nm,'mRNA_perArea','maxproj_dir','spots_dir','prob_thresh','pxl_to_micron','receivers_label','chan_nm');

disp('done')


%%
mat_dir = 'D:\2021-01-01-smFISH\Dapi_Lefty_Produce_Nodal\MIP\WT_jux_WT_wnti';%WT_jux_NKO
close all
ff = dir(mat_dir);
cd(mat_dir);
all_dat = [];
q = 1;
dat = [];
receive = [];
produce = [];
d_mean=struct;
d_err =struct;
p_mean=[];
all_dat=[];
for jj=1:size(ff,1)
    if ~isdir(ff(jj).name) &&  ~isempty(regexp(ff(jj).name,'WT_jux_WT_wnti_' ,'ONCE')) && ~isempty(strfind(ff(jj).name,'.mat'))
        disp(ff(jj).name)
        load(ff(jj).name);
        all_dat = [];
        for jj=1:size(mRNA_perArea,2)
            all_dat=cat(2,mRNA_perArea.rcv);
            p_mean = mean(all_dat(:,2:2:end),2);
            p_err = std(all_dat(:,2:2:end),[],2);
            vect1 = all_dat(:,1);
        end
        
        receive(q).dat = p_mean;
        receive(q).err  = p_err;
         figure(2), errorbar(vect1,receive(q).dat,receive(q).err,'-p','LineWidth',2); hold on%,dat(:,2)

        q = q+1;
        
    end
   
end

legend('20 hr','4 hr','9 hr');
ylabel('N of Nodal mRNA spots per um^2 in receivers')
xlabel('Distance from producing cells, um')
title('WT producers -WT receivers, in presence of WNTi')
ylim([0 0.008])

%%

close all
for jj=1:size(mRNA_perArea,2) 
figure(1),plot(mRNA_perArea(jj).prdc(:,1),mRNA_perArea(jj).prdc(:,2),'b');hold on
end
xlabel('Distance from producing cells, um')
ylabel('Number of mRNA spots per um^2')
h = figure(1);
h.CurrentAxes.FontSize = 14;
h.CurrentAxes.LineWidth = 2;
%ylim([0 max(expr_dat_inreceiving(:,1))+100]);
title('Nodal mRNA spots')
legend('in producing cells')
for jj=1:size(mRNA_perArea,2) 
figure(2),plot(mRNA_perArea(jj).rcv(:,1),mRNA_perArea(jj).rcv(:,2),'r');hold on
end
title('Nodal mRNA spots')
legend('in receiving cells')


