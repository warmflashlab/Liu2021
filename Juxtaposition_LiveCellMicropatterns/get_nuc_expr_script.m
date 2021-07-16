%% get expression in nyclei of cells, no two-cell types
close all
clear all
direc = 'D:\2020-11-07-nodal\plate2_Dapi_blank_blank_pSmad2\E6_SB';%D:\2020-09-02-DapiNodalmembraneLefty_40X\40x\20200628_003_20200628_95055 PM\max_proj
cd(direc)
ff = readAndorDirectory(direc);
chan_nm = {'DAPI','pSmad2'}; %order as in chan_tmp 
% channel order:  [dapi ; stain_to_quantify; producing_cells_marker ] indexing starts from 1
chan_tmp = [1 4];
thresh2 = 0.9;
smallstuff = 250;
usedapi = 1;
norm_chan=1;%indexed into chan_tmp
raw_img=struct;
pxl_to_micron = 0.325;%20X=0.65, 40x = 0.325 on SD confocal. on LSM : 20X = 0.759, 40X = 
img_fin = struct;
save_nm= [ direc 'pSmad2_E6_SB.mat'];
chan = 2;
positions = struct;
for i=1:size(ff.p,2)
pos=ff.p(i);%i-1
fnm=getAndorFileName(ff,pos,[],[],ff.w(chan_tmp(1)));
if pos<10
ilastik2 = [direc '\plate2_f000' num2str(pos) '_MIP_w0000_Probabilities.h5'];
[masks]=readIlastikProbMask(ilastik2,thresh2);%
else
ilastik2 = [direc '\plate2_f00' num2str(pos) '_MIP_w0000_Probabilities.h5']; 
[masks]=readIlastikProbMask(ilastik2,thresh2);%
end

img_fin = get_maxpr_raw_img(chan_tmp,ff,pos,masks);
[expr_dat_all]=get_nuc_expr(img_fin,masks,chan,norm_chan);
 positions(i).nodist = expr_dat_all;
end

save(save_nm,'positions','chan_nm','chan_tmp','direc','expr_dat_all','masks')
%% plot

close all
clear all
expr_thresh_high=10; % remove data points that have normalized expression higher than this thresh
fn = 'D:\2020-11-07-nodal\plate2_Dapi_blank_blank_pSmad2\E6pSmad2_E6.mat';
[e6_cells,e6_err]=get_plot_dat(fn,expr_thresh_high);

fn = 'D:\2020-11-07-nodal\plate2_Dapi_blank_blank_pSmad2\E6_SBpSmad2_E6_SB.mat';
[e6_sb,e6_sb_err]=get_plot_dat(fn,expr_thresh_high);
dat = [];
dat(:,1) = cat(2,e6_cells,e6_sb);
dat(:,2) = cat(2,e6_err,e6_sb_err);
figure(1), errorbar(1:2,dat(:,1),dat(:,2),'p','Markersize',15,'Linewidth',2,'Color','m');hold on
title('pSMAD2 in  cell nuclei')
ylabel('pSmad2/DAPI')
xlim([0 3])
ylim([0.2 0.6])
h = figure(1);
xticks([1,2]);
xticklabels({'E6 ','E6 + SB'});
h.CurrentAxes.FontSize = 14;
h.CurrentAxes.LineWidth = 2;
xtickangle(-20)

