function [status] = write_modem_covariance(cfile,cov,xsmooth,ysmooth,zsmooth,nsmooth)
%
%  Usage:  [status] = writeCov_3D(cfile,cov,xsmooth,ysmooth,zsmooth,nsmooth)
%  
%  Write a 3-D covariance matrix & simple rules in the Mod3DMT format.
%  Optional smoothing arguments define the model covariance.
%
% Inputs:
%       -cfile is a filename string
%       -cov is an array of covariance flag integers
%       -xsmooth is a [nz x 1] vector of smoothing in the x direction for each layer
%       -ysmooth is a [nz x 1] vector of smoothing in the y direction for each layer
%       -zsmooth is a single number between 0 and 1 for vertical smoothing
%       -nsmooth is the number of times smoothing is applied
%
%  In the covariance matrix, air = 0 and ocean = 9.
%
%  (c) Anna Kelbert, 2009

[nx,ny,nz] = size(cov);
fid = fopen(cfile,'w');


%  initialize smoothing parameters
if nargin < 3
    xsmooth = 0.3 * ones(1,nz);
    ysmooth = 0.3 * ones(1,nz);
    zsmooth = 0.3;
end

%  how many times the smoothing should be applied 
if nargin < 6
    nsmooth = 1;
end

%  write the fixed header
header = [ ...
'+-----------------------------------------------------------------------------+';...
'| This file defines model covariance for a recursive autoregression scheme.   |';...
'| The model space may be divided into distinct areas using integer masks.     |';...
'| Mask 0 is reserved for air; mask 9 is reserved for ocean. Smoothing between |';...
'| air, ocean and the rest of the model is turned off automatically. You can   |';...
'| also define exceptions to override smoothing between any two model areas.   |';...
'| To turn off smoothing set it to zero. This header is 16 lines long.         |';...
'| 1. Grid dimensions excluding air layers (Nx, Ny, NzEarth)                   |';...
'| 2. Smoothing in the X direction (NzEarth real values)                       |';...
'| 3. Smoothing in the Y direction (NzEarth real values)                       |';...
'| 4. Vertical smoothing (1 real value)                                        |';...
'| 5. Number of times the smoothing should be applied (1 integer >= 0)         |';...
'| 6. Number of exceptions (1 integer >= 0)                                    |';...
'| 7. Exceptions in the form e.g. 2 3 0. (to turn off smoothing between 2 & 3) |';...
'| 8. Two integer layer indices and Nx x Ny block of masks, repeated as needed.|';...
'+-----------------------------------------------------------------------------+';...
];
for i = 1:16
    fprintf(fid,'%s\n',header(i,:));
end
fprintf(fid, '\n');

% record # 1: size of the model
fprintf(fid,'%d %d %d\n',nx,ny,nz);
fprintf(fid, '\n');

% record # 2: X-smoothing
count = 0;
for k = 1:nz
    fprintf(fid,'%G ',xsmooth(k));
    count = count+1;
    if rem(count,10)==0
        fprintf(fid, '\n');
    end
end
fprintf(fid, '\n');
fprintf(fid, '\n');

% record # 3: Y-smoothing
count = 0;
for k = 1:nz
    fprintf(fid,'%G ',ysmooth(k));
    count = count+1;
    if rem(count,10)==0
        fprintf(fid, '\n');
    end
end
fprintf(fid, '\n');
fprintf(fid, '\n');

% record # 4: Z-smoothing
fprintf(fid,'%G\n',zsmooth);
fprintf(fid, '\n');

% record # 5: number of times to apply smoothing operator
fprintf(fid,'%d\n',nsmooth);
fprintf(fid, '\n');

% records # 6 & 7: edit exceptions manually if needed
fprintf(fid,'%d\n',0);
fprintf(fid, '\n');

% record # 8: covariance matrix in the convention "get what you see"
fprintf(fid, '\n');
for k = 1:nz
    if rem(k,10)==0 || k==1 || k==nz
        disp(['Layer ',num2str(k),' ... ...'])
        if k==nz
            disp('Finished saving covariance file')
        end
    end
   fprintf(fid,'%d %d\n',k,k);
   for i = nx:-1:1
       % y index incremented fastest
       for j = 1:ny
        fprintf(fid,'%d ',cov(i,j,k));
       end    
       fprintf(fid, '\n');
    end
end

% close file
status = fclose(fid);

