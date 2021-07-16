% get LSM Nodal movies

%read the laser scanning data and get max projections
% save max projections
img_subfolder=struct;
    q = 1;  
    tp = 0;

for xx =3% 2:16;
    direc = ['C:\Users\Asya\Dropbox (Warmflash Lab)\LabFiles\Cecilia\200908_YuxtaposNodalLefty\Track000' num2str(xx) '\'];
    direc2save = 'D:\2020-09-08-Cecilia_JuxtaposNodalLefty\';
    
    ff = dir(direc);
    % get the number of time points and the names of the subfoilders to
    % open
    for jj=1:size(ff,1)
        if ~isempty(regexp(ff(jj).name,'.files','ONCE'))
            tp = tp+1;
            img_subfolder(q).nm = ff(jj).name(1:length(ff(jj).name)-6);
            q = q+1;
        end
    end 
    for kk=1:tp
        
    for chan = 1:2       
        multitp_nuc = [];        
        reader = bfGetReader(img_subfolder(kk).nm);
        multitp_nuc = bfMaxIntensity(reader,1,chan);
        %figure(1), imshow(multitp_nuc,[]);
        disp(['populating  tp ' num2str(kk)]);
            if (xx-1)<10
                imwrite(multitp_nuc,[direc2save '\' 'Nodal_p000' num2str(xx-1) '_w000' num2str(chan) '.tif'],'writemode','append');
            end
            if (xx-1)>=10
                imwrite(multitp_nuc,[direc2save '\' 'Nodal_p00' num2str(xx-1) '_w000' num2str(chan) '.tif'],'writemode','append');
            end
        
        disp('done')
    end
    end
end
    