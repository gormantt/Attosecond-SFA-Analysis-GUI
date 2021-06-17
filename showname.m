clear;
filelist =dir('*.sif');
[m n]=size(filelist);
s_f=zeros(m,1);
for ii = 1:m
    s_f(ii)= length(filelist(ii).name);
end
filenames=zeros(m,max(s_f),'uint8');
filenames = char(filenames);
for ii = 1:m
    filenames(ii,1:s_f(ii))= filelist(ii).name;
end
if exist('namelist.txt','file')
    delete('namelist.txt');
end
fid = fopen('namelist.txt','w');
fprintf(fid,'%d\n',m);
fprintf(fid,'%d\n',max(s_f));
for ii=1:m
    fprintf(fid,'%s\n',filenames(ii,:));
%     fprintf(fid,'\n');
end
fclose(fid);