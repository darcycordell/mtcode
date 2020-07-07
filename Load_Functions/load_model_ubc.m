function [m,d] = load_model_ubc(name,meshname,mesh_corner)
%
% Function to load a UBC model into the standard model format
%
% Usage: [m,d] = load_model_ubc(name,meshname,mesh_corner)
%
% Inputs:
%   name: the name of the model file. This file should contain one column
%   of values (e.g. resistivity or density)
%
%   meshname: the name of the mesh file. This file should contain 5 lines
%       Line 1: ny nx nz
%       Line 2: UTM origin of the SW top corner of the mesh with elevation
%       in meters above sea level (+ numbers above sea level)
%       Line 3: The y thicknesses (e.g. 69*250 means 69 cells with 250 thickness)
%       Line 4: The x thicknesses (e.g. 43*100 42*150 means 43 cells with
%               100 m thickness followed by 42 cells with 150 m thickness
%       Line 5: The z thicknesses (same formats as above)
%
%   mesh_corner: You need to manually convert the UTM origin (line 2 of the
%   mesh file) into latitude and longitude. This can be done here:
%       http://www.rcn.montana.edu/Resources/Converter.aspx
%
%   For example, suppose the origin in the UBC format .mesh file is 
%   Easting = 355000 and Northing = 5999000 in UTM zone 19S.
%   These are converted to lat long with:
%       mesh_corner = [-36.128465 -70.61138200735813]
%
% Outputs:
%       "m" is a standard model structure
%       "d" is a standard data structure


%
% For debugging purposes:
% name = 'SimPEG_inv_l0l2_3.0_0.2_0.025.den';
% meshname = 'ubc_mesh_expanding.mesh';

mod = load(name);
m.name = name;
m.niter = '';


fid = fopen(meshname); if fid == -1; error('Mesh file not found'); end;

tline = 1; i=1;
while tline ~= -1

    tline = fgetl(fid);
    line{i} = tline;
    i = i+1;
end

m.ny = str2double(line{1}(1:2));
m.nx = str2double(line{1}(4:5));
m.nz = str2double(line{1}(7:8));

utm = strsplit(line{2},' ');

m.origin(2) = 0; 
m.origin(1) = 0; 
m.origin(3) = -str2double(utm{3});

%Read in blocks from mesh file
y = strsplit(line{3},' ');
x = strsplit(line{4},' ');
z = strsplit(line{5},' ');

%Get cell thicknesses from blocks
m.dx = [];
for i = 1:length(x)
    xx = strsplit(x{i},'*');
    m.dx = [m.dx; str2double(xx{2})*ones(str2double(xx{1}),1)];
end

m.dy = [];
for i = 1:length(y)
    yy = strsplit(y{i},'*');
    m.dy = [m.dy; str2double(yy{2})*ones(str2double(yy{1}),1)];
end

m.dz = [];
for i = 1:length(z)
    zz = strsplit(z{i},'*');
    m.dz = [m.dz; str2double(zz{2})*ones(str2double(zz{1}),1)];
end

%Vector of cell edge locations starting from zero
m.x = [0; cumsum(m.dx)]+m.origin(1);
m.y = [0; cumsum(m.dy)]+m.origin(2);
m.z = [0; cumsum(m.dz)]+m.origin(3);

%Find cell center midpoints
m.cx = (m.x(1:end-1)+m.x(2:end))/2;
m.cy = (m.y(1:end-1)+m.y(2:end))/2;
m.cz = (m.z(1:end-1)+m.z(2:end))/2;

%Create a meshgrid of cell centers
[m.X,m.Y]=meshgrid(m.cy,m.cx);

m.A = reshape(mod,[m.nz m.nx m.ny]);

m.A = permute(m.A, [3 2 1]);

m.A(m.A==-100) = NaN;

m.npad = [0 0];
%%
[m.lon,m.lat] = utm2geo(m.cy+500000,m.cx,mesh_corner(2),mesh_corner(1));

%Build meshgrid of longitudes and latitudes.
[m.LON,m.LAT]=meshgrid(m.lon,m.lat);  
% 
% plot(m.LON(:),m.LAT(:),'.k'); hold on; plot_geoboundaries
%%
%The centered Easting = 363500 (east_center_utm) and centered Northing =
%6007500 (north_center_utm) in UTM zone 19S.
d = make_nan_data;

d.origin(1) = (max(m.lat)+min(m.lat))/2;
d.origin(2) = (max(m.lon)+min(m.lon))/2;


d.origin(2) = mesh_corner(2);
d.origin(1) = mesh_corner(1);

[m.y,m.x] = geo2utm(m.lon,m.lat,d.origin(2),d.origin(1));
m.y = m.y-500000;



%Determine elevation topography surface
m.Z = zeros(m.nx,m.ny);
for i = 1:m.nx
    for j = 1:m.ny
        
        ind = find(isnan(squeeze(m.A(i,j,:))),1,'last');
        
        if isempty(ind)
            ind = 0;
        end
        
        m.Z(i,j) = m.cz(ind+1);
        
    end
end

m.origin(2) = m.y(1); 
m.origin(1) = m.x(1); 


end

