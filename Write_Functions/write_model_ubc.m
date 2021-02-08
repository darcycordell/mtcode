function write_model_ubc(modfile,m)
%
% Function which writes UBC format mesh and conductivity files.
%
% Usage: write_model_ubc(modfile,m)write_model_ubc(modfile,m)
%
% Inputs:
%    modfile is the name of the model file to write
%    m is Matlab structure containing model information

% get cell indices
indx = 1:m.nx;
indy = 1:m.ny;
indz = 1:m.nz;

% get cell widths
dx = m.dx;
dy = m.dy;
dz = m.dz;

% get number of cells
nx = m.nx;
ny = m.ny;
nz = m.nz;

% get south west corner of mesh
x_south = min(m.x);
y_west = min(m.y);
z_top = min(m.z); % min since m.z is in metres b.s.l

%%
%Mesh file cell order is EW, NS, Z, therefore need to swap x and y from above
% also need to negate m.z
disp('Writing mesh file')
fid2=fopen([modfile,'.mesh'], 'w');
fprintf(fid2, '%d %d %d \n', ny, nx, nz);
fprintf(fid2, '%d %d %d \n', round(y_west), round(x_south), -round(z_top));

for i=1:ny
    fprintf(fid2, '%d ', round(dy(i)));
end
fprintf(fid2, '\n');
for i=1:nx
    fprintf(fid2, '%d ', round(dx(i)));
end
fprintf(fid2, '\n');
for i=1:nz
    fprintf(fid2, '%d ', round(dz(i)));
end
fclose(fid2);

%%
disp('Writing .con file')
m.A(isnan(m.A))=10^17;
A = 1./m.A; % convert to conductivity (S/m)

fid3=fopen([modfile,'.con'], 'w');
% according to UBC documentation, depth index changes fastest, then
% easting, then northing, i.e. in order from innermost loop to outer loop:
for i=indx % then move to next cell to north
    for j=indy % then output cells to east
        for k=indz % first output every cell in column
            fprintf(fid3, '%E\n ', A(i,j,k));
        end
    end
end
fclose all;

end