%% lefty nodal system:
close all
clear all
% 1. get masks of the secreting cells from dapi+mamebrane
% 2. get masks of the receiver cells ( as not producer cells)
direc='D:\2021-01-01-Lefty_2-4-9-20hr\40x\4_9_20h\20hr\LKO_LKO_Act_20hr';%D:\2020-09-02-DapiNodalmembraneLefty_40X\40x\20200628_003_20200628_95055 PM\max_proj
cd(direc)
sv_dir = 'D:\2021-01-01-Lefty_2-4-9-20hr\40x\4_9_20h\20hr\';
ff = readAndorDirectory(direc);
chan_nm = {'DAPI','membrane','Membrane'}; %order as in chan_tmp 
save_nm= [ sv_dir 'mmb_LKO_LKO_Act_20hr.mat'];
manual_bg=[];%600 720 620
% channel order:  [dapi ; stain_to_quantify; producing_cells_marker ] indexing starts from 1
tosave=0;
chan_tmp = [1 3 3];
receivers_mmb=0;        % NOTE which cells have the membrane marker, receivers or producers!!!!!!!!!!!!!!!!!
z = [];                 % set this to specific z plane if want to look at Lefty in specific image plane
direc_z = '';           %D:\2020-09-02-DapiNodalmembraneLefty_40X\40x\20200628_003_20200628_95055 PM
thresh2 = 0.92;
smallstuff =6500;
smallstuff2=4000;%
usedapi =0;%if true, use dapi in combination with membrane marker to isolate producing cells, and then receiving cells
norm_chan=[];% if norm channel is not empty, then also need to have Dapi masks (quantification is pxl-based but for nuclei only)

raw_img=struct;
pxl_to_micron = 0.325;%  20X=0.65, 40x = 0.325 on SD confocal. on LSM : 20X = 0.759, 40X = 
q_dist = 2; % channel for which to quantify the distance dependence (indexed into chan_tmp )
dxy = 10;   % in microns, 'slices' to look at away from the cell border  
d_max= 300; % in microns, how far from border to look
img_fin = struct;
stats_tmp=[];
inreceiving=1; % if set to 1, then the quantification is in the receiving cells with respect to producing cells (if 0, then swapped)
% input params: 
% get_expression_pxl(images,pxl_to_micron,producing_cells_mask,dxy,d_max,receiving_cells_mask,chan_to_quantify,normalization
% channel,smallstuff)
positions = struct;
%v = [1 2 3 6 9];
for i=1:size(ff.p,2)%size(v,2)%%%%%
pos=ff.p(i);%v(i)-1;%%;%;%;%;%;%i-1
if receivers_mmb
[mask_produce_fin,mask_receiving] = get_masks(pos,direc,chan_tmp,thresh2,smallstuff,usedapi);
else
[mask_receiving,mask_produce_fin] = get_masks(pos,direc,chan_tmp,thresh2,smallstuff,usedapi);
end

%-----only used for images, where there are borders with empty space
% stats_tmp = regionprops(mask_produce_fin,'Area');
% mask_produce_fin=bwareaopen(mask_produce_fin,max(cat(1,stats_tmp.Area)));
%------------
imshowpair(mask_receiving,mask_produce_fin)

if inreceiving
    if usedapi
        img_fin = get_maxpr_raw_img(chan_tmp,ff,pos,mask_receiving | mask_produce_fin,manual_bg);
    else
        img_fin = get_maxpr_raw_img(chan_tmp,ff,pos,mask_receiving,manual_bg);
    end
    [expr_dat,dist_i,expr_dat_all]=get_expression_pxl(img_fin,pxl_to_micron,mask_produce_fin,dxy,d_max,mask_receiving,q_dist,norm_chan,smallstuff2);
else
    if usedapi
        img_fin = get_maxpr_raw_img(chan_tmp,ff,pos,mask_receiving | mask_produce_fin,manual_bg);
    else
        img_fin = get_maxpr_raw_img(chan_tmp,ff,pos,mask_produce_fin,manual_bg);
    end
    [expr_dat,dist_i,expr_dat_all]=get_expression_pxl(img_fin,pxl_to_micron,mask_receiving,dxy,d_max,mask_produce_fin,q_dist,norm_chan,smallstuff2);
end
 figure(10), errorbar(dist_i(2:end),expr_dat(:,1),expr_dat(:,2));hold on
 positions(i).dat = expr_dat;
 positions(i).nodist = expr_dat_all;
end
dist_vect = dist_i; % in microns
if inreceiving 
title([ chan_nm{q_dist} ' in receiving cells'])
if ~isempty(norm_chan)
    title([ chan_nm{q_dist} ' /DAPI in receiving nuclei'])
end
xlabel('Distance from Nodal-producing cells')
else
    title([ chan_nm{q_dist} ' in producing cells'])
    if ~isempty(norm_chan)
    title([ chan_nm{q_dist} ' /DAPI in receiving nuclei'])
end
    xlabel('Distance from receiving cells')
end
ylabel('Fluorescence intensity') 
if tosave
save(save_nm,'positions','dist_vect','dxy','d_max','mask_produce_fin','mask_receiving','pxl_to_micron','chan_nm','chan_tmp','direc','inreceiving','q_dist','expr_dat_all')
end
figure,imshowpair(img_fin(3).dat,mask_receiving) %mask_produce_fin pixels that are considered receiving cells/areas
figure, imshow(img_fin(3).dat,stretchlim(img_fin(3).dat))%
close all
all_dat = cat(2,positions.dat);
% get aver over positions

p_mean_ko = mean(all_dat(:,1:2:end),2);
p_err_ko = std(all_dat(:,1:2:end),[],2);
figure(2), errorbar(dist_vect(2:end),p_mean_ko,p_err_ko,'b');hold on

%% plot just the average expression, no distnace dependence
close all
%clear all
expr_thresh_high=10; % remove data points that have normalized expression higher than this thresh
fn = 'D:\2020-11-07-nodal\max_proj\bcatcells_DapiMembr_Smad2\pSmad2_inproducing.mat';
[wt_produce,wt_err_produce]=get_plot_dat(fn,expr_thresh_high);


%% plot averaged over positions as a function of distance from border 
close all
hold on
load('D:\2020-12-23-lefty\40x\MIP\matfiles\Lefty_inreceiving_LKO_6hr.mat')%NODALmRNA_inreceiving_KO    MIPLeftymRNA_inreceiving_KO
all_dat = cat(2,positions.dat);
% get aver over positions

p_mean_ko = mean(all_dat(:,1:2:end),2);
p_err_ko = std(all_dat(:,1:2:end),[],2);
load('D:\2020-12-23-lefty\40x\MIP\matfiles\Lefty_inreceiving_LKO_ACT_6hr.mat')%MIPLeftymRNA_inreceiving_wt NODALmRNA_inreceiving_WT   
all_dat = cat(2,positions.dat);
p_mean_wt = mean(all_dat(:,1:2:end),2);
p_err_wt = std(all_dat(:,1:2:end),[],2);

figure(2), errorbar(dist_vect(2:end),p_mean_ko,p_err_ko,'b');hold on
figure(2), errorbar(dist_vect(2:end),p_mean_wt,p_err_wt,'r');hold on
title('Expression in receiving cells, 20X (4 images)')
if inreceiving 
title(['Mean ' chan_nm{q_dist} ' in receiving cells'])
xlabel('Distance from Nodal-producing cells')
else
    title([ chan_nm{q_dist} ' in producing cells'])
    xlabel('Distance from receiving cells')
end
ylabel('Fluorescence intensity')
h = figure(2);
h.CurrentAxes.FontSize = 14;
h.CurrentAxes.LineWidth = 2;
legend('LKO-ACT','LKO + ACT');

%% plot mean expression vs time for time course experiments
clear all
mat_dir = 'D:\2021-01-01-smFISH\Dapi_Lefty_Produce_Nodal\MIP\20hr_lefty_mRNA\';%D:\2021-01-01-smFISH\Dapi_Lefty_Produce_Nodal\MIP\WT_jux_WT
close all
clc
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
    if ~isdir(ff(jj).name) &&  ~isempty(regexp(ff(jj).name,'_WT_jux_WT' ,'ONCE')) && ~isempty(strfind(ff(jj).name,'.mat'))%meanLefty_WT_WT_
        disp(ff(jj).name)
        load(ff(jj).name);
        all_dat = cat(2,positions.dat);
        d_mean(q).dat = mean(all_dat(:,1:2:end),2,'omitnan');
        d_err(q).dat = std(all_dat(:,1:2:end),[],2);
        receive(q,1) = mean(cat(1,positions.nodist),'omitnan');
        receive(q,2)  = std(cat(1,positions.nodist),'omitnan');%/power(size(cat(1,positions.nodist),1),0.5);
        q = q+1;        
    end
end
vect1 = [1:size(receive,1)];

 figure(2), errorbar(vect1,receive(:,1),receive(:,2),'pr','LineWidth',2); hold on%,dat(:,2)
% ylabel('Mean pixel intensity (LEFTY)')
% xlabel('Time point, hours')
% title('20X; Average over positions')
% h = figure(2);
% h.CurrentAxes.LineWidth = 2;
% h.CurrentAxes.FontSize = 12;
%ylim([0 350])
dist_dat = cat(2,d_mean.dat);
dist_err = cat(2,d_err.dat)./power(size(dist_dat,2),0.5);%
for jj=1:size(vect1,2)
figure(3), errorbar(dist_vect(2:end)',dist_dat(:,jj),dist_err(:,jj),'-p','LineWidth',2); hold on%,dat(:,2)
end
legend('WT-WT','WT-WTnoind','WT-WT-WNTi');%'WT-WT','WT-WTnoind','WT-WT-WNTi','NKO-WT','WT-NKO','WT-NKO-noind'
ylabel('Lefty mRNA (mean pixel intensity)')
xlabel('Distance from producing cells, um')
title('Mean Lefty mRNA, 20 hr time point ')
h = figure(3);
h.CurrentAxes.LineWidth = 2;
h.CurrentAxes.FontSize = 12;
ylim([0 300])

