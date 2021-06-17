function [signal back ref] = SifSBR(filename)
fid=fopen(filename);
%Confirm this is an appropriate Andor file.
FileType=fgetl(fid);
if ~isequal(FileType,'Andor Technology Multi-Channel File')
   fclose(fid);
   error('Not the expected Andor SIF image file.');
end
%Verify appropriate version
FileVersion=readword(fid);
if ~logicor(str2double(FileVersion),65538)
   fclose(fid);
   error('Unexpected file version. May work, but needs verification.');
end
[variables fileinfo]=SifDetails(fid);
signal = SifFrame(fid,fileinfo,1);
[back backend]= SifBack(fid,fileinfo);
ref = SifRef(fid,fileinfo,backend);
fclose(fid);
