function [back backend]= SifBack(fid,fileinfo)
Npixels=fileinfo.dimX*fileinfo.dimY;

%Adjust the position to the beginning of the desired frame. 
newpos=fileinfo.datastart+Npixels*(fileinfo.Nframes)*4;

%Adjust the position to the desired positions
fseek(fid,newpos,'bof');
tlin=fgetl(fid);
if tlin==('0')
    tlin=fgetl(fid);
    if tlin~=('1')
        fclose(fid);
        error('end of file');
    end
elseif tlin~=('1');
    fclose(fid);
    error('end of file');
end     

for ii = 1:22
    tlin=fgetl(fid);
    if ~ischar(tlin)
        fclose(fid);
        error('Wrong file');
    end
end
back=fread(fid,Npixels,'float32=>float32');

%Convert it into a matrix, and then transpose it to orient it normally. 
back=reshape(back,fileinfo.dimX,fileinfo.dimY)';
back=double(back);
backend=ftell(fid);