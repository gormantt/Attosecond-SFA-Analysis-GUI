function ref = SifRef(fid,fileinfo,backend)
Npixels=fileinfo.dimX*fileinfo.dimY;

%Adjust the position to the desired positions
fseek(fid,backend,'bof');
for ii = 1:23
    fgetl(fid);
end
ref=fread(fid,Npixels,'float32=>float32');

%Convert it into a matrix, and then transpose it to orient it normally. 
ref =reshape(ref,fileinfo.dimX,fileinfo.dimY)';
ref=double(ref);
