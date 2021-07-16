function [expr_dat_all]=get_nuc_expr(img_fin,masks,chan,norm_chan)

mask = imfill(masks,'holes');% prediff mask
expr_dat_all=[];
stats_all = regionprops(mask,img_fin(chan).dat,'Centroid','MeanIntensity');

if norm_chan
stats_norm = regionprops(mask,img_fin(norm_chan).dat,'Centroid','MeanIntensity');    
expr_dat_all = cat(1,stats_all.MeanIntensity)./cat(1,stats_norm.MeanIntensity);
else
expr_dat_all = cat(1,stats_all.MeanIntensity);    
end