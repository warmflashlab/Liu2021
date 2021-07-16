% process the movie from the Nodal-Lefty imaging
% no dapistain, only membrane producing cells stain
% todo: 
% 1. load the movie images and the ilasitk prob mask
% 2. make the masks for receiving and producing cell areas
% 3. quantify Nodal in Lefty cells and in Nodal cells in time

maxproj_dir='D:\2020-10-13-inserts_Noggin_SBmovies\normal';
ff = dir(maxproj_dir);
segm_chan = 2;% channel used to get the segementation (masks of receiving/producing cells areas)
cd(maxproj_dir);
mask_producing = struct;
mask_receiving= struct;
fnm = struct;
ilastik = struct;
dt =15; % time interval between movie frames
prob_thresh=0.9;
small_stuff=500;
positions = [2];%1 2
quantify_chan = 1; % which channel to use for getting expression data (for  Nodal =  chan 1), indedex as w0001 : chan   = 1
dyn_dat = struct;
% params for kymograph for nodal movie
pxl_to_micron = ;%20X ( pxltomicron = 0.325 - 40X)
t_i=1;
t_end=2;
dxy = 2;% size of slices to look at incrementally, away from the borderin microns
max_d = 50;%  max distance from the border ,in microns
dat = [];
kymo_dat=struct;
in_nodal_cells =0; % if set to 1, get data for producing cells
for pos =3%1:size(positions,2)    
    q = 1;
for jj=1:size(ff,1)
    if ~isdir(ff(jj).name) &&  ~isempty(regexp(ff(jj).name,['p000' num2str(pos-1) '_w000' num2str(segm_chan) ],'ONCE')) && ~isempty(strfind(ff(jj).name,'.tif'))
        fnm(q).name = ff(jj).name;% this is the image/movie name
        ilastik(q).name = [fnm(q).name(1:end-4) '_Probabilities.h5'];
        disp(ilastik(q).name);
        mask_producing.alltimes = imfill(bwareaopen(readIlastikProbMask(ilastik(q).name,prob_thresh),small_stuff),'holes');
        mask_receiving.alltimes =  bwareaopen(~mask_producing.alltimes,small_stuff);
        [expression_dat,nT,img_dat]=nodal_dyn(pos-1,quantify_chan,ff,mask_receiving, mask_producing);
        dyn_dat(pos).receiving = cat(1,expression_dat.receive);
        dyn_dat(pos).producing = cat(1,expression_dat.produce);     
        %------------------- 
       [spacetime]=get_kymo(t_i,t_end,img_dat,pxl_to_micron,mask_receiving,mask_producing,in_nodal_cells,dxy,max_d) ;
        %----------------------------
        
        q = q+1;
        kymo_dat(pos).dat = spacetime;

    end  
    
end

end
%save('nodal_movie_dat_inproducingcells','kymo_dat','dyn_dat','maxproj_dir','pxl_to_micron','small_stuff','quantify_chan','in_nodal_cells','segm_chan');
disp('done');
%% plot nodal dynamics in each cell type
%load('nodal_movie_dat_inproducingcells.mat');
t_end = 79;
t_vect = (1:t_end )*dt/60;
close all
receive_all=[];
produce_all = [];
for jj=1:size(dyn_dat,2)
    receive_all = cat(2,receive_all,dyn_dat(jj).receiving(1:t_end ));
    figure(3), plot(dyn_dat(jj).receiving(1:t_end ));hold on;
    produce_all = cat(2,produce_all,dyn_dat(jj).producing(1:t_end ));
    figure(4), plot(dyn_dat(jj).producing(1:t_end ));hold on;

end
r = mean(receive_all,2);
err_r = std(receive_all,[],2);
p = mean(produce_all,2);
err_p = std(produce_all,[],2);

figure(1), errorbar(t_vect(1:2:end),r(1:2:end),err_r(1:2:end),'LineWidth',2);hold on
errorbar(t_vect(1:2:end),p(1:2:end),err_p(1:2:end),'LineWidth',2);
legend('Nodal in receiving cells','Nodal in producing cells');
ylabel('Nodal intensity')
xlabel('Time, hrs')
h1 = figure(1);
h1.CurrentAxes.LineWidth = 2;
h1.CurrentAxes.FontSize = 14;
title('Nodal LSM movie, 20X, 3 positions')

%% plot kymographs

load('nodal_movie_dat_inproducingcells.mat');
load('nodal_movie_dat.mat');

close all
mean_kymo = kymo_dat(1).dat ;%((+ kymo_dat(2).dat+kymo_dat(3).dat))/3
        figure(10), s = pcolor(mean_kymo);hold on
        xlabel('Time, frames')
        if in_nodal_cells
            ylabel('Distance from Receiving cells, microns')
            title('Nodal in producing cells')
        else
            ylabel('Distance from producing cells, microns')
            title(' Nodal in receiving cells')
        end
        colorbar
        h  = figure(10);
        s.FaceColor = 'interp';
        set(s,'EdgeColor','none');
        
       % caxis([0 1500]);


