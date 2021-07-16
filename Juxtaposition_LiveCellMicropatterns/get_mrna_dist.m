function [expr_dat,dist_i]=get_mrna_dist(pxl_to_micron,mask_producing,dxy,d_max,img_mrna)
 expr_dat=[];
mask = mask_producing;% prediff mask

I_tmp = bwdist(mask);
%imshow(I_tmp);
dist_img_um = I_tmp*pxl_to_micron; % convert image distance to microns (now each pixel is how many microns away is the nearest prediff pixel)
dist_i = (1:dxy:d_max);% in microns
expr_dat=[];
q = 0;
for w=1:size(dist_i,2)-1
pxls_bin = find(dist_img_um >=(dist_i(w)+q)& dist_img_um <=dist_i(w+1));% linear index of pixels
[(dist_i(w)+q) dist_i(w+1)];
% now look at the values of these pixels in the receiver image
test_img = zeros(size(mask));
test_img(pxls_bin) = 1;
nuc_pxl = [];
nuc_pxl = (test_img & img_mrna);
stats_n = [];
stats_n = regionprops(nuc_pxl,'Centroid');
stats_area = regionprops(test_img,'Area');
tmp_dat = [];    
    if ~isempty(stats_n)
        tmp_dat = size(stats_n,1)/(stats_area.Area/(pxl_to_micron*pxl_to_micron)); % number of  mRNA spots in that slice away from border, normalized to area
        if w == 0
            figure(1), imshowpair(imdilate(nuc_pxl,strel('disk',1)),mask);hold on%masks mask_pluri
        end
    end    
  if ~isempty(tmp_dat)
  tmp_dat(tmp_dat==Inf)=nan;
  expr_dat(w,1) = tmp_dat;
  end
  
end

end