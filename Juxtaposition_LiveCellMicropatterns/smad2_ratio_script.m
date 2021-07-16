%% get nuclear to cyto ratio in static images, lefty-nodal system, no mask for cyto
% 1. use maxprojsctipr to generate indexed and regular max projections
close all
clear all
pos1 =16;%13-16
pluri_chan = [0 1];
maxproj_dir = 'D:\2021-01-01-Jux_4_9_20hr\Smad2_plate2_20210101_70456 PM\MIP\';
indx_img_dir='D:\2021-01-01-Jux_4_9_20hr\Smad2_plate2_20210101_70456 PM\MIPidx\';
chan_tmp = [1 2 3];% channel order:  [dapi ; stain_to_quantify; producing_cells_marker ] indexing starts from 1
chan_nm = {'Dapi','Smad2/3','MembraneProducers'};
clear masks
thresh2 = 0.92;
smallstuff =400;
inreceivers = 1;
receiver_mmb=0;% which cell type is labelled with membrane
if pos1 <10
prefix_str = 'Smad2_plate2_MIP_f000';%'normal-SB-noggin-nc-c-TLC-S4-ESC_MIP_p000'
else
prefix_str = 'Smad2_plate2_MIP_f00';    
end
usedapi=1;
[mask_produce_fin,mask_receiving] = get_masks(pos1-1,maxproj_dir,chan_tmp,thresh2,smallstuff,usedapi);
%figure, imshowpair(mask_produce_fin,mask_receiving);
masks = struct;
if receiver_mmb
masks(1).alltimes = imdilate(mask_produce_fin,strel('disk',1));
masks(2).alltimes = imdilate(mask_receiving,strel('disk',1));
figure,imshowpair(mask_receiving,imfill(masks(2).alltimes,'holes'));
else
%mask_receiving = mask_produce_fin;    
masks(1).alltimes = imdilate(mask_receiving,strel('disk',1));
masks(2).alltimes = imdilate(mask_produce_fin,strel('disk',1));
figure,imshow(imfill(masks(2).alltimes,'holes'));
end
% todo: need to get better masks for separate nuclei
%% get indexed nucler image
pos = pos1-1;
cd(indx_img_dir)
%indx_nm = ['p000' num2str(pos) 'MIPidx.tif'];
ff = dir(indx_img_dir);
for jj=1:size(ff,1)
    if pos <10
if ~isdir(ff(jj).name) && ~isempty(strfind(ff(jj).name,['Smad2_plate2_f000' num2str(pos) ]))     
    indx_nm = ff(jj).name;    
end
    else
       if ~isdir(ff(jj).name) && ~isempty(strfind(ff(jj).name,['Smad2_plate2_f00' num2str(pos) ]))    
    indx_nm = ff(jj).name;    
       end
    end        
end
reader = [];
indx_im=struct;
    reader = bfGetReader(indx_nm);% reader for specific channel
    nz=reader.getSizeZ;  
    nT = reader.getSizeT; 
    chan = 1;
    for t = 1
     iPlane=reader.getIndex(nz - 1, chan -1, t - 1) + 1;
        indx_im(t).dat=bfGetPlane(reader,iPlane);
        figure(1),imshow(indx_im(t).dat,[]);
    end
idx_img=indx_im(t).dat ;    
%% get raw z-stack
%tic
spacetime=[];
pxl_to_micron = 0.325;
dxy = 5;%15
cyto_sz=6;
dt = [];    % min between frames
d_fin = 200; % microns, how many microns to do overall (from interface)
str = ['Smad2_plate2_f00' num2str(pos) '_w000' num2str(chan_tmp(2)-1)];
z_slice_dir ='D:\2021-01-01-Jux_4_9_20hr\Smad2_plate2_20210101_70456 PM';% [maxproj_dir '\all_chans']; % path to the directory with all the time groups and zslices 

[img]=get_z_planes(z_slice_dir,str);
bg = [620];%bg for the quantified channel
%toc
t = [];
[ratio,dist_i]=get_ratio_pxlrange(idx_img,pxl_to_micron,masks,dxy,cyto_sz,img,t,bg,d_fin);
% d = cat(2,dist_i,dist_i(end)+dxy )';
 dat = cat(2,dist_i(2:end)',ratio(:,1));
figure, plot(dat(:,1),dat(:,2));
save([indx_img_dir '\ratio_dat_NKO-WT-noind-' num2str(pos) '.mat'],'dat','cyto_sz','dxy','inreceivers','chan_tmp','bg','chan_nm','pxl_to_micron','ratio','dist_i','pos','maxproj_dir');
disp('done')    
    
%% plot results
clc
close all
mat_dir = 'D:\2021-01-01-Jux_4_9_20hr\Smad2_plate1_20210101_70441 PM\MIPidx';
str = 'ratio_dat_WT-WT_p';
[dist_dat,dist_err,dist_vect,str_nm]=get_mean_dist(mat_dir,str);
figure(3), errorbar(dist_vect,dist_dat,dist_err,'-p','LineWidth',2); hold on%,dat(:,2)
ylabel('Nuclear:cyto ratio')
xlabel('Distance from producing cells, um')
title(['WT-WT: 20 hr time point, N:C ratio of ' str_nm 'in receivers']);
h = figure(3);
h.CurrentAxes.LineWidth = 2;
h.CurrentAxes.FontSize = 11;   
    