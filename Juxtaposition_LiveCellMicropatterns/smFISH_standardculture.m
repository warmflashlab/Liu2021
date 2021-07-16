close all
maxproj_dir='D:\2020-12-01-smFISH\StandardCulture_FISH_dapiNodalLefty\smFISH_quantification\';
ff = dir(maxproj_dir);
segm_chan = 1;% dapi
q_chan =3;% [dapi nodal lefty] [1 2 3]
chan_nm = {'Dapi','Nodal','Lefty'};
cd(maxproj_dir);
mask=[];
fnm = struct;
ilastik = struct;
prob_thresh=0.9;
% params for kymograph for nodal movie
pxl_to_micron = 0.3250;%( pxltomicron = 0.325 - 40X)
small_stuff=500;
mean_expr=[];
v = [1 2 4 8];
dat = [];
for tp=1:size(v,2)   
for pos =1%:size(positions,2)    

        fnm.name =  ['Activin' num2str(v(tp)) 'h_MIP_p000' num2str(pos-1) '_w000' num2str(segm_chan-1) '.tif' ];
        img = imread(['Activin' num2str(v(tp)) 'h_MIP_p000' num2str(pos-1) '_w000' num2str(q_chan-1) '.tif' ]);
        img_no_bg = img-mean(img(~imdilate(mask,strel('disk',30))));
        ilastik.name = [fnm.name(1:end-4) '_Probabilities.h5'];
        disp(ilastik.name);
        mask = imfill(bwareaopen(readIlastikProbMask(ilastik.name,prob_thresh),small_stuff),'holes');                 
        imshowpair(img,mask);
stats = regionprops(imdilate(mask,strel('disk',10)),img_no_bg,'Area','Centroid','MeanIntensity');
%imshow(mask);hold on
%xy = cat(1,stats.Centroid);
%text(xy(:,1),xy(:,2),num2str(cat(1,stats.Area)));
mean_expr = mean(cat(1,stats.MeanIntensity));
% load the .csv file manually (for now)
spots_data = importdata([maxproj_dir num2str(v(tp)) '_p000' num2str(pos-1) '_w000' num2str(q_chan-1) '.csv']);
disp(['Loaded : ' num2str(v(tp)) '_p000' num2str(pos-1) '_w000' num2str(q_chan-1) '.csv,  cnannel: ' chan_nm{q_chan}]);
spots_data(spots_data==0)=NaN;
img_mrna = zeros(size(mask));
for jj=1:size(spots_data,1)
    if ~isnan(spots_data(jj,:))
img_mrna(spots_data(jj,1),spots_data(jj,2))=1;
    end
end
%imshowpair(img_mrna,mask)
%hold on
%plot(spots_data(:,2),spots_data(:,1),'p','LineWidth',0.05);hold on
%[~]=get_mrna_percell(mask,img_mrna); % still to write this

stats_n = [];
stats_n = regionprops(mask,'Area'); % how many cells are there on the image
avg_cell_area = mean(cat(1,stats_n.Area));
cells = sum(cat(1,stats_n.Area))/avg_cell_area;
n_mRNA = size(spots_data,1); % how many total mRNA spots were detected
RNA_per_cell = n_mRNA/cells;

end
dat(tp,1)=v(tp);
dat(tp,2) = RNA_per_cell ;
dat(tp,3) = cells;
dat(tp,4) = mean_expr;

end
%close all
figure(6), plot(dat(:,1),dat(:,2),'-rp','LineWidth',2);hold on
xlabel('Time, hrs')
ylabel('mRNA molecules per cell')
title(chan_nm{q_chan})
h6 = figure(6);
h6.CurrentAxes.FontSize = 12;
h6.CurrentAxes.LineWidth = 2;


ylim([min(dat(:,2))-10 max(dat(:,2))+10])
figure(7), plot(dat(:,1),dat(:,3),'-gp','LineWidth',2);hold on
xlabel('Time, hrs')
ylabel('Cells per image')
title(chan_nm{q_chan})
ylim([min(dat(:,3))-10 max(dat(:,3))+10])
h7 = figure(7);
h7.CurrentAxes.FontSize = 12;
h7.CurrentAxes.LineWidth = 2;



