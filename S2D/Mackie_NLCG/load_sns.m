function [ny,nz,dy,dz,sens]=load_sns(name,version)
%==========================================================================
fid=fopen(name,'r');

A=fscanf(fid,'%i',2);
ny=A(1);
nz=A(2);
dy=fscanf(fid,'%f',ny);
dz=fscanf(fid,'%f',nz);
if version==6
   fscanf(fid,'%i',1);
end
sens=fscanf(fid,'%f',ny*nz);
sens=reshape(sens,ny,nz)';
sens=log10(sens);

fclose(fid);
end