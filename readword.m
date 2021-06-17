function out=readword(fid)

%Read each word, removing trailing white space
tmp=textscan(fid,'%s ',1);
out=tmp{:};