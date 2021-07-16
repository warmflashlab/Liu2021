function [dist_dat,dist_err,dist_vect,str_nm]=get_mean_dist(mat_dir,str)
ff = dir(mat_dir);
cd(mat_dir);
q = 1;
d_mean=struct;
for jj=1:size(ff,1)
    if ~isdir(ff(jj).name) &&  ~isempty(regexp(ff(jj).name,str ,'ONCE')) && ~isempty(strfind(ff(jj).name,'.mat'))
        disp(ff(jj).name)
        load(ff(jj).name);
        d_mean(q).dat = dat(:,2);            
        q = q+1;        
    end
end    
dist_vect = dat(:,1);
dist_dat = mean(cat(2,d_mean.dat),2);
dist_err = std(cat(2,d_mean.dat),[],2);
str_nm = chan_nm{chan_tmp(2)};

end