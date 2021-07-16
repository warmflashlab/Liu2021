function [pxl_dat,dist_i]=sign_range(img_dat,pxl_to_micron,mask_receiving,mask_producing,time,dxy,max_d)
 
% get the pixel intensity depending on distance from the border between
% cells
% for nodal movies
pxl_dat = [];
mask = mask_producing.alltimes(:,:,time);% producing cells
I_tmp = bwdist(mask);
dist_img_um = I_tmp*pxl_to_micron; % convert image distance to mucrons (now each pixel is how many microns away is the nearest prediff pixel)
%figure, imshow(dist_img_um,[]);
dist_i = (1:dxy:max_d);% in microns
pxl_dat=[];
q = 0;
for w=1:size(dist_i,2)-1
pxls_bin = find(dist_img_um >=(dist_i(w)+q)& dist_img_um <=dist_i(w+1));% linear index of pixels
[(dist_i(w)+q) dist_i(w+1)];
% now look at the values of these pixels in the receiving cells
test_img = zeros(size(mask));
test_img(pxls_bin) = 1;
nuc_pxl = [];
nuc_pxl = (test_img & mask_receiving.alltimes(:,:,time));
%figure(1), imshowpair(nuc_pxl,img_dat(time).nobg);hold on
stats_n = [];
stats_n = regionprops(nuc_pxl,img_dat(time).nobg,'Centroid','MeanIntensity','PixelIdxList');
pxl_dat(w) = mean(nonzeros(cat(1,stats_n.MeanIntensity)));

end
end
