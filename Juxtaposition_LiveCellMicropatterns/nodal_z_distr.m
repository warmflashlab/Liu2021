%% extract Nodal distribution along z-planes

close all
clear all
% TODO: needs to be normalized per something or taken per cell...
direc_z = 'D:\2020-08-24-DapiNodalMembrane_Lefty';%D:\2020-09-02-DapiNodalmembraneLefty_40X\40x\20200628_003_20200628_95055 PM\max_proj
direc_il = 'D:\2020-08-24-DapiNodalMembrane_Lefty\max_pr';
cd(direc_il)
ff_z = readAndorDirectory(direc_z);
ff_il = readAndorDirectory(direc_il);
chan_nm = {'DAPI','Lefty','Membrane'}; %order as in chan_tmp 
% channel order:  [dapi ; stain_to_quantify; producing_cells_marker ] indexing starts from 1
chan_tmp = [1 4 3];
z = []; % set this to specific z plane if want to look at Lefty in specific image plane
thresh2 = 0.9;
smallstuff =200;
pxl_to_um=0.617284;
um_per_plane=0.9231;%0.9091;
nz = 11;
ratio_fin=struct;
producing = 0;
for i=1:size(ff_il.p,2)
pos=ff_il.p(i);%i-1
fnm_dapi  =  getAndorFileName(ff_il,pos,[],[],ff_il.w(chan_tmp(1)));
ilastik_dapi = [fnm_dapi(1:(length(fnm_dapi)-4)) '_Probabilities.h5'];
fnm_mmb  =  getAndorFileName(ff_il,pos,[],[],ff_il.w(chan_tmp(3)));
ilastik_mmb = [fnm_mmb(1:(length(fnm_mmb)-4)) '_Probabilities.h5'];

[mask_all]=readIlastikProbMask(ilastik_dapi,thresh2);%
[mask_mmb]=readIlastikProbMask(ilastik_mmb,thresh2);%
mask_nuc = imfill(bwareaopen(mask_all,smallstuff),'holes');
mask_mmb = imfill(mask_mmb,'holes');
mask_mmb_2 = bwareaopen(mask_mmb&~mask_nuc,smallstuff);

tmp = imfill(imdilate(mask_mmb_2,strel('disk',20)),'holes');
producing_nuc=tmp&mask_nuc;
if producing
imshowpair(producing_nuc,mask_mmb_2);
end
if producing == 0
receiving = mask_all&~producing_nuc;
receiving = bwareaopen(receiving,smallstuff);
receiving =imfill(receiving,'holes');
dil=imdilate(receiving,strel('disk',5));
mask_mmb_2 = bwareaopen(dil&~receiving,smallstuff);  
 
 producing_nuc=receiving;
 imshowpair(producing_nuc,mask_mmb_2);
end
%close all,figure,imshowpair(receiving,mask_mmb_2)
expr_mmb=zeros(nz,1);
expr_nuc=zeros(nz,1);
ratio=zeros(nz,1);
for z=1:nz
[img_fin,nz]=get_z_raw_img(chan_tmp,ff_z,pos,z,mask_mmb_2| mask_all);
stats_mmb = regionprops(mask_mmb_2,img_fin(2).dat,'Area','MeanIntensity');
stats_nuc = regionprops(producing_nuc,img_fin(2).dat,'Area','MeanIntensity');
if z==5
    figure(1), imshow(img_fin(2).dat,[]);
end
expr_mmb(z,:) = mean(cat(1,stats_mmb.MeanIntensity));%
expr_nuc(z,:) = mean(cat(1,stats_nuc.MeanIntensity));
ratio(z,:)=expr_mmb(z,:)/expr_nuc(z,:);
end
ratio_fin(i).dat = ratio;
end

x_vect = (1:nz)*um_per_plane;
m_ratio = mean(cat(2,ratio_fin.dat),2);
err = std(cat(2,ratio_fin.dat),[],2);
figure(1),errorbar(x_vect,m_ratio,err,'Color','r','LineWidth',2);hold on
%plot(expr_nuc);
%legend('membrane','nuc')
xlabel('Microns (Basal->Apical)')
if producing
title([ chan_nm{2} '  localization in producing cells']);
else
title([ chan_nm{2} '  localization in receiving cells']);    
end
ylabel(['Membrane-to-Nuclear ratio of '  chan_nm{2}]);
h = figure(1);
h.CurrentAxes.FontSize = 18;
h.CurrentAxes.LineWidth = 3;

for jj=1:size(ratio_fin,2)
figure(1), plot(ratio_fin(jj).dat); hold on
end


