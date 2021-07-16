function [img_fin]=get_maxpr_raw_img(chan_tmp,ff,pos,mask_allnuc,manual_bg)
 raw_img=struct;
 img_fin=struct;
 bg=[];
for jj=1:size(chan_tmp,2)
    if ispc
        fntmp = getAndorFileName(ff,pos,[],[],ff.w(chan_tmp(jj)));
           % fntmp = [ num2str(pos) '_MIP_w000' num2str(chan_tmp(jj)-1) '.tif'];
       
        raw_img(jj).dat =imread(fntmp);%%%%%%% jj
        disp(fntmp);
    else
        raw_img(jj).dat =imread(fntmp);%%%%%%% jj
        disp(fntmp);
        
    end
end
for jj=1:size(chan_tmp,2) 
       if isempty(manual_bg)
        bg(jj) = mean(raw_img(jj).dat(mask_allnuc));% mean value of pixels in the area of receiving cells
       else
        bg(jj) =  manual_bg(jj);
       end
        disp(bg(jj));
        
        %img_fin(jj).dat = simplebg([],mask_allnuc,raw_img(jj).dat);%chan_tmp(jj)
        img_fin(jj).dat = raw_img(jj).dat-bg(jj);%chan_tmp(jj)
        figure(jj), imshow(img_fin(jj).dat,[]);%img_fin(jj).dat        
end
