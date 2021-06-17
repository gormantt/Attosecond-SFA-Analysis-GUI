fid = fopen('namelist.txt','r');
m=str2double(fgetl(fid));
n=str2double(fgetl(fid));
filenames=zeros(m,n,'uint8');
filenames = char(filenames);
for ii=1:m
    filenames(ii,:) = fgetl(fid);
end
fclose(fid);