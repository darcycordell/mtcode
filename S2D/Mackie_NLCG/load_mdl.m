function [ny,nz,dy,dz,rho]=load_mdl(name,version)
%==========================================================================

disp('IN : load_mdl')

fid=fopen(name,'r');

A=fscanf(fid,'%i',2);
ny=A(1);
nz=A(2);
dy=fscanf(fid,'%f',ny);
dz=fscanf(fid,'%f',nz);
if version==6
   fscanf(fid,'%i',1);
end
rho=fscanf(fid,'%f',ny*nz);
rho=reshape(rho,ny,nz)';
rho=log10(rho);

fclose(fid);
end