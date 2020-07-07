function C = load_covariance(filename,model_filename)
%Function which loads a given ModEM covariance file and places it within a
%standard model structure format where A is the covariance mask integers.
%
% Inputs:
%   filename - The ModEM covariance filename
%   model_filename - The corresponding model filename associated with the
%      covariance
%
% Outputs:
%   C - A standard model structure which contains all the info such that
%   the covariance masks can be plotted using standard functions such as
%   plot_slice
%       The C structure also contains additional info about the covariance
%       parameters stored in the cov_params variable.

m = load_model_modem(model_filename);

fid = fopen(filename);

line = fgetl(fid);

%Read header
while line(1) == '+' || line(1) == '#' || line(1) == '|'
    line = fgetl(fid);
    if isempty(line)
        line = fgetl(fid);
    end
end

dim = str2num(line);
nx = dim(1); ny = dim(2); nz = dim(3);

line = fgetl(fid);

cov_params.covx = fscanf(fid,'%f',nz);
cov_params.covy = fscanf(fid,'%f',nz);
cov_params.covz = fscanf(fid,'%f',1);

cov_params.num_times_smoothing_applied = fscanf(fid,'%i',1);
cov_params.num_exceptions = fscanf(fid,'%i',1);

while ~strcmp(line,'1 1')
    line = fgetl(fid);
end

A = zeros(nx,ny,nz);
for i = 1:nz
   block = fscanf(fid,'%i',nx*ny);
   block = flipud(block);
   A(:,:,i) = reshape(block,ny,nx)';
   A(:,:,i) = A(:,end:-1:1,i);
   line = fscanf(fid,'%i',2);
end

fclose(fid);


A(A==0) = NaN;

if nx==m.nx && ny==m.ny && nz==m.nz && isequal(size(A), size(m.A)) || (isvector(A) && isvector(m.A) && numel(A) == numel(m.A))
    C = m;
    C.A = A;
    C.cov_params = cov_params;
else
    error('The size of the covariance model space is not the same as your input model space. Check file inputs');
end
