function [site,locxy,ref] = load_sitefile(name,ns)
%Function which loads the site names from the *.stn file from a 2D NLCG
%Inversion
%
% Usage: [site,locxy,ref] = load_sitefile(name,ns)

fid=fopen(name,'r');

fgetl(fid);
ref = str2double(fgetl(fid));

i=1;
site{ns} = []; locxy = ones(ns,3);

while 1
    
    line = fgetl(fid);
    
    if line~=-1
    site{i} = line(1:find(line==';')-1);
    
    locxy(i,2) = str2double(line(find(line==';')+1:end)); % *1000; think it's already in m
    
    i=i+1;
    else
        break
    end


end

end