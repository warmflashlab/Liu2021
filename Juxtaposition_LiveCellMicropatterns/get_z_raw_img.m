function [img_fin,nz]=get_z_raw_img(chan_tmp,ff,pos,z,mask_allnuc)
 raw_img=struct;
 img_fin=struct;
 bg = [];
for jj=1:size(chan_tmp,2)
    if ispc
        fntmp = getAndorFileName(ff,pos,[],[],ff.w(chan_tmp(jj)));
        [img_z,nz]=get_z_plane(fntmp,z);
        raw_img(jj).dat =img_z;
       % disp(fntmp);
    else
        [img_z,nz]=get_z_plane(fntmp,z);
        raw_img(jj).dat =img_z;
        disp(fntmp);
        
    end
end
for jj=1:size(chan_tmp,2)     
        bg(jj) = mean(raw_img(jj).dat(~mask_allnuc));% get values of pixels that  don't have nuclei
        %disp(bg(jj));
        %img_fin(jj).dat = simplebg([],mask_allnuc,raw_img(jj).dat);%chan_tmp(jj)
        img_fin(jj).dat = raw_img(jj).dat-bg(jj);%chan_tmp(jj)
        %figure(jj), imshow(img_fin(jj).dat,[]);%img_fin(jj).dat        
end
