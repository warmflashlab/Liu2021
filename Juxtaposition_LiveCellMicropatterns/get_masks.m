function [mask_receiving,mask_produce_fin] = get_masks(pos,direc,chan_tmp,thresh2,smallstuff,usedapi)
%
mask_produce_fin=[];
ff = readAndorDirectory(direc);
if usedapi % if true, use dapi in combination with membrane marker to isolate producing cells, and then receiving cells
fnm_dapi  =   getAndorFileName(ff,pos,[],[],ff.w(chan_tmp(1)))% ilastik file name for dapi channel  
ilastik_dapi = [fnm_dapi(1:(length(fnm_dapi)-4)) '_Probabilities.h5'];
[mask_dapi]=readIlastikProbMask(ilastik_dapi,thresh2);% 
mask_dapi = imfill(mask_dapi,'holes');
fnm_produce = getAndorFileName(ff,pos,[],[],ff.w(chan_tmp(3)-1));% chan_tmp(3)-1) only tmp, for when the MIP contains only segmented data
ilastik2 = [fnm_produce(1:(length(fnm_produce)-4)) '_Probabilities.h5'];
[mask_nodal]=readIlastikProbMask(ilastik2,thresh2);% mask nodal producing cells
mask_nodal = imfill(mask_nodal,'holes');
imshowpair(mask_dapi,imfill(mask_nodal,'holes'));

mask_produce_fin = mask_dapi & mask_nodal;
mask_produce_fin = bwareaopen(mask_produce_fin,smallstuff);
% final
mask_receiving =  bwareaopen(mask_dapi & ~mask_produce_fin,smallstuff);
%figure, imshowpair(mask_receiving,mask_produce_fin);


else
fnm_produce = getAndorFileName(ff,pos,[],[],ff.w(chan_tmp(3)));% ilastik file name for the producing cells
ilastik2 = [fnm_produce(1:(length(fnm_produce)-4)) '_Probabilities.h5'];
ilastik2
[mask_nodal]=readIlastikProbMask(ilastik2,thresh2);% mask nodal producing cells
mask_nodal = imfill(imdilate(mask_nodal,strel('disk',10)),'holes');
%imshowpair(mask_allnuc,imfill(mask_nodal,'holes'))
%figure, imshow(imfill(mask_allnuc,'holes') & imfill(mask_nodal,'holes'))
mask_produce_fin = bwareafilt(mask_nodal,[smallstuff 1000000]);
% final
mask_receiving =  bwareaopen(~mask_produce_fin,smallstuff);
%figure, imshowpair(mask_receiving,mask_produce_fin);

end
end