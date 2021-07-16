function smoothMask = makeColonyMaskUsingDapiImage(image1, metadata)
%% ---------------- returns colonymask -----------------------------
%% -----------------------------------------------------------------
%%
% convert dapi to binary
minI = double(min(image1(:)));
maxI = double(max(image1(:)));
forIlim = mat2gray(image1);

t = 0.8*graythresh(forIlim)*maxI + minI;
debris = image1 > 8*t; % super bright non-cell stuff
debris = imdilate(debris, strel('disk', 15));

iFactor = 1; %[ for images with high signal to background (too bright), keep iFactor >=1, else keep iFactor from 0-1]
dapi2 = image1 > iFactor*t;
dapi2 = dapi2&~debris;
%figure; imshow(dapi2);

% keep only the colony
xres = double(metadata.xres);
s = round(20/xres); % assuming a cell diameter is 10 microns
dapi3 = imclose(dapi2, strel('disk', round(s/1.2)));
dapi3 = imfill(dapi3,'holes');
%dapi3 = imopen(dapi3, strel('disk', round(s/3)));
dapi3 = imerode(dapi3, strel('disk',round(s/2))); % play around with the radius
cleanMask = bwareafilt(dapi3, 1, 'largest');
%figure; imshow(cleanMask);

% smooth edges
windowSize = 200;
kernel = ones(windowSize) / windowSize ^ 2;
blurryMask = conv2(single(cleanMask), kernel, 'same');
smoothMask = blurryMask > 0.5;
smoothMask = imdilate(smoothMask, strel('disk', round(s/2)));
%figure; imshowpair(image1, smoothMask);

end