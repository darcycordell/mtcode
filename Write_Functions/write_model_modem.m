function write_model_modem(outputfile,dx,dy,dz,A,origin)
%
% Function which writes a ModEM model file format
%
% Usage: write_model_modem(outputfile, dx, dy, dz, A, origin)
%
% Inputs:
%   -outputfile is a string filename
%   -dx is a vector cell thicknesses in the north-south direction
%   -dy is a vector of cell thicknesses in the east-west direction
%   -dz is a vector of cell thicknesses in the vertical direction
%   -A is a (nx x ny x nz) array of resistivity values
%   -origin is a 3 x 1 vector specifying the southwest top corner of the
%   model in model coordinates [south west top]

%DC July 2015: Function to save output model in ModEM format which can be
%read into ModEM inversion

A(isnan(A)==1)=10^17;
A=log(A);
type='LOGE';

nzAir = 0; %DC: I am not sure how ModEM handles air layers. My understanding is that there are zero air layers
rotation=0;

write_WS3d_model(outputfile,dx,dy,dz,A,nzAir,type,origin,rotation);

end %END MAIN


function status = write_WS3d_model(fname,x,y,z,rho,nzAir,type,origin,rotation,niter)
% writes a 3D resistivity model in Weerachai Siripunvaraporn's format;
% allows for natural log resistivity by setting type = 'LOGE'
%  (c) Anna Kelbert, 2009
%  open file for output
fid = fopen(fname,'w');
[nx, ny, nz] = size(rho);
if nargin <= 6
    type = ' ';
end
% output file
%comment = ['Written by Matlab write_WS3d_model script Iteration No.',num2str(niter+1,'% 4.0f')]; 
if nargin < 10
    fprintf(fid, '# Written by Matlab write_WS3d_model script\n');
else
    fprintf(fid, '# Written by Matlab write_WS3d_model script Iteration No.% 4.0f\n', niter+1);
end
fprintf(fid, '%d %d %d %d %s\n', nx, ny, nz, 0, type);
for j = 1:nx
    status = fprintf(fid,'%G ',x(j));
end
fprintf(fid, '\n');
for j = 1:ny
    status = fprintf(fid,'%G ',y(j));
end
fprintf(fid, '\n');
for j = 1:nz
    status = fprintf(fid,'%G ',z(j));
end
fprintf(fid, '\n');
for k = 1:nz
    if rem(k,10)==0 || k==1 || k==nz
        disp(['Layer ',num2str(k),' ... ...'])
        if k==nz
            disp('Finished saving model file')
        end
    end
   fprintf(fid,'\n');
   for j = 1:ny
       % x index incremented fastest
       for i = nx:-1:1
        fprintf(fid,'%15.5E',rho(i,j,k));
       end    
       fprintf(fid, '\n');
    end
end
%  add origin, rotation angle to end
if nargin <= 7
    origin = [-sum(x)/2 -sum(y)/2 0];
    rotation = 0;
end
if ~isnan(origin)
    fprintf(fid, '%d %d %d\n', origin);
    fprintf(fid, '%d\n', rotation);
end
status = fclose(fid);
end