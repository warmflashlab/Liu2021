function [expr_dat,dist_i]=get_mrna_percell(mask,img_mrna)
stats_n = [];
stats_n = regionprops(mask,'PixelIdxList','Centroid'); % how many cells are there on the image
stats_area = regionprops(test_img,'Area'); %
tmp_dat = [];   






end

