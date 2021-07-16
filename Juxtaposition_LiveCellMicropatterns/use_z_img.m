function [img_fin]=use_z_img(chan_tmp,ff,pos,mask_allnuc)
 raw_img=struct;
 img_fin=struct;
 bg=[];
for jj=1:size(chan_tmp,2)
    if ispc
        fntmp = getAndorFileName(ff,pos,[],[],ff.w(chan_tmp(jj)));
        raw_img(jj).dat =imread(fntmp);%%%%%%% jj
        disp(fntmp);
    else
        raw_img(jj).dat =imread(fntmp);%%%%%%% jj
        disp(fntmp);
        
    end
end
for jj=1:size(chan_tmp,2) 
       
        bg(jj) = mean(raw_img(jj).dat(~mask_allnuc));% get values iof pixels that  that don'e have nuclei
        disp(bg(jj));
        %img_fin(jj).dat = simplebg([],mask_allnuc,raw_img(jj).dat);%chan_tmp(jj)
        img_fin(jj).dat = raw_img(jj).dat-bg(jj);%chan_tmp(jj)
        figure(jj), imshow(img_fin(jj).dat,[]);%img_fin(jj).dat        
end
