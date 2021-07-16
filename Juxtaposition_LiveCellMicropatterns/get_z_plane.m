function [img_z,nz]=get_z_plane(fnm,z)
img_z=[];
reader = bfGetReader(fnm);
nz=reader.getSizeZ;
chan = 1;
time = [];
for jj = z
            %iPlane=reader.getIndex(jj - 1, chan -1, time - 1) + 1;
            img_z=bfGetPlane(reader,z);
           % figure(jj),imshow(img_z,[]);            
            
end
end