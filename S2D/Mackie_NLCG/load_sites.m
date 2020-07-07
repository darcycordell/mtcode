function [indsit,sitrms]=load_sites(name)
%==========================================================================
% get the number of sites and index

fid=fopen(name,'r');
line=1; 
m=0;
while line ~= -1
    line=fgetl(fid);
    k=findstr(line,'block');
    if ~isempty(k)
        m=m+1;
        indsit(m)=sscanf(line(k+5:max(size(line))),'%i');
        line=fgetl(fid); line=fgetl(fid);
        word=strread(line,'%s');
        sitrms(m)=str2num(char(word(7)));
    end
end
fclose(fid);
return