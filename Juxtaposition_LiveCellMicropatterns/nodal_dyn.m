function [expression_dat,nT,img_dat]=nodal_dyn(pos,quantify_chan,ff,mask_receiving, mask_producing,pos_str,overlap)
q = 1;
fnm2=struct;
img_dat = struct;
expression_dat=struct;
bg= [];
for jj=1:size(ff,1)
    if ~isdir(ff(jj).name) &&  ~isempty(regexp(ff(jj).name,[pos_str num2str(pos) '_w000' num2str(quantify_chan) ],'ONCE')) && ~isempty(strfind(ff(jj).name,'.tif'))
        fnm2(q).name = ff(jj).name;% this is the image/movie name
        disp(['Movie used for expression quantification: ' num2str(ff(jj).name)])
        reader = bfGetReader(fnm2(q).name); % fname needs to be changed,
        nT = reader.getSizeT;
        nz = 1; % since using max projections
        chan =1;
        % make a structure with all image time points
        img_z=[];
        tmp_bg_mask=[];
        for h=1:nT-1
            iPlane=reader.getIndex(nz - 1, chan -1, h-1) + 1;
            img_z=bfGetPlane(reader,iPlane);
            bg(h) = mean(img_z(mask_receiving.alltimes(:,:,1))); % mean pixel intensity in the area of receiving cells at time point 1 (sicne no signaling is happening)
            img_dat(h).nobg = img_z-bg(h);
            
            if overlap
                
              im = img_dat(h).nobg;%
              targetsize = size(im)-overlap;
              r = centerCropWindow2d(size(im),targetsize);
              img_dat(h).nobg = imcrop(im,r);  
                
            end
            
            
             figure(2), imshowpair(mask_receiving.alltimes(:,:,h),mask_producing.alltimes(:,:,h));hold on%
            %figure(1), imshow(img_dat(h).nobg,[]);  hold on
            % this is Nodal in producing cells, no distance dependence
            tmp1 = [];
            tmp1 = regionprops(mask_producing.alltimes(:,:,h),img_dat(h).nobg,'MeanIntensity','Area');
            tmp2 = [];
            tmp2 = regionprops(mask_receiving.alltimes(:,:,h),img_dat(h).nobg,'MeanIntensity','Area');
            expression_dat(h).produce = mean(nonzeros(cat(1,tmp1.MeanIntensity)));
            expression_dat(h).receive = mean(nonzeros(cat(1,tmp2.MeanIntensity)));
        end
        q = q+1;% this is incrementing the positions
        
    end
end
