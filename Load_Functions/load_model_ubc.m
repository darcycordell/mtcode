function [m,d] = load_model_ubc(name,meshname,mesh_corner)
%
% Function to load a UBC model into the standard model format
%
% This function assumes the model contains contains values of 
% electrical conductivity. In order to convert to the standard m structure, 
% this function converts to resistivity and sets all values greater than 
% 1e15 to NaN.
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
%       http://rcn.montana.edu/Resources/Converter.aspx
%
%   For example, suppose the origin in the UBC format .mesh file is 
%   Easting = 355000 and Northing = 5999000 in UTM zone 19S.
%   These are converted to lat long with:
%       mesh_corner = [-36.128465 -70.61138200735813]
%
% Outputs:
%       "m" is a standard model structure
%       "d" is a standard data structure, but only containing the variable
%       d.origin (center of mesh)


%
% For debugging purposes:
% name = 'SimPEG_inv_l0l2_3.0_0.2_0.025.den';
% meshname = 'ubc_mesh_expanding.mesh';

mod = load(name); % assumed conductivity file?
m.name = name;
m.niter = '';

fid = fopen(meshname); if fid == -1; error('Mesh file not found'); end

tline = 1; i=1;
while tline ~= -1
    tline = fgetl(fid);
    line{i} = tline;
    i = i+1;
end

yxz = str2num(line{1});
m.ny = yxz(1);
m.nx = yxz(2);
m.nz = yxz(3);

utm = strsplit(line{2},{' ','\t'});
if length(utm)==5
    utm(1)=[];
    utm(end) = [];
end
m.origin(2) = str2double(utm{1}); 
m.origin(1) = str2double(utm{2}); 
m.origin(3) = -str2double(utm{3}); % negate for m b.s.l. convention

%Read in blocks from mesh file
y = strsplit(line{3},' ');
x = strsplit(line{4},' ');
z = strsplit(line{5},' ');

if y(end)==""
    disp('Extra space ignored at end of line 3')
    y(end)=[];
end
if x(end)==""
    disp('Extra space ignored at end of line 4')
    x(end)=[];
end
if z(end)==""
    disp('Extra space ignored at end of line 5')
    z(end)=[];
end

fclose(fid);

% Get cell thicknesses from blocks
% check if file is written in format nx*dx below
m.dx = [];
for i = 1:length(x)
    if isempty(regexp(x{i},'*', 'once'))
        m.dx = [m.dx; str2double(x{i})];
    else
        xx = strsplit(x{i},'*');
        m.dx = [m.dx; str2double(xx{2})*ones(str2double(xx{1}),1)];
    end
end

m.dy = [];
for i = 1:length(y)
    if isempty(regexp(y{i},'*', 'once'))
        m.dy = [m.dy; str2double(y{i})];
    else
        yy = strsplit(y{i},'*');
        m.dy = [m.dy; str2double(yy{2})*ones(str2double(yy{1}),1)];
    end
end

m.dz = [];
for i = 1:length(z)
    if isempty(regexp(z{i},'*', 'once'))
        m.dz = [m.dz; str2double(z{i})];
    else
        zz = strsplit(z{i},'*');
        m.dz = [m.dz; str2double(zz{2})*ones(str2double(zz{1}),1)];
    end
end

%Vector of cell edge locations starting from zero
m.x = [0; cumsum(m.dx)]+m.origin(1);
m.y = [0; cumsum(m.dy)]+m.origin(2);
m.z = [0; cumsum(m.dz)]+m.origin(3);

%Find cell center midpoints
m.cx = (m.x(1:end-1)+m.x(2:end))/2;
m.cy = (m.y(1:end-1)+m.y(2:end))/2;
m.cz = (m.z(1:end-1)+m.z(2:end))/2;

% create meshgrid of cell boundaries
[m.X,m.Y]=meshgrid(m.y,m.x);
%Create a meshgrid of cell centers
[m.Xc,m.Yc]=meshgrid(m.cy,m.cx);

% m.A = reshape(1./mod,[m.nz m.nx m.ny]);
% m.A = permute(m.A, [3 2 1]);

% slower method below but clearer what/where values are being assigned
mod = 1./mod; % convert from conductivity to resistivity
m.A = zeros(m.nx,m.ny,m.nz);
count = 0;
for ix = 1:m.nx
    for iy = 1:m.ny
        for iz = 1:m.nz
            count = count + 1;
            m.A(ix,iy,iz) = mod(count);
        end
    end
end

% m.A(m.A==-100) = NaN;
m.A(m.A>1e15) = NaN; % assuming these values represent air cells! might need different condition for other physical proprties

m.npad = [sum(m.dx~=min(m.dx))/2 sum(m.dy~=min(m.dy))/2]; % assuming same number of positive and negative paddings cells
%%
if m.ny ~= m.nx
    [m.lon,~] = utm2geo(m.cy+500000,m.cx(1),mesh_corner(2),mesh_corner(1));
    [~,m.lat] = utm2geo(m.cy(1)+500000,m.cx,mesh_corner(2),mesh_corner(1));
else
    [m.lon,m.lat] = utm2geo(m.cy+500000,m.cx,mesh_corner(2),mesh_corner(1));
end

%Build meshgrid of longitudes and latitudes.
[m.LON,m.LAT]=meshgrid(m.lon,m.lat);  

%%
% debugging: The centered Easting = 363500 (east_center_utm) and centered Northing =
%6007500 (north_center_utm) in UTM zone 19S.
d = make_nan_data;

d.origin(1) = (max(m.lat)+min(m.lat))/2;
d.origin(2) = (max(m.lon)+min(m.lon))/2;
% 
% debugging: these are calculated above... not needed anymore?
% d.origin(2) = mesh_corner(2);
% d.origin(1) = mesh_corner(1);
% [m.y,m.x] = geo2utm(m.lon,m.lat,d.origin(2),d.origin(1));
% m.y = m.y-500000;

%Determine elevation topography surface
m.Z = zeros(m.nx,m.ny);
for i = 1:m.nx
    for j = 1:m.ny
        
        ind = find(isnan(squeeze(m.A(i,j,:))),1,'last');
        
        if isempty(ind)
            ind = 0;
        end
        
        if ind == m.nz
            ind = m.nz-1;
        end
        
        m.Z(i,j) = m.z(ind+1);
        
    end
end

% m.origin(2) = m.y(1); % not needed?
% m.origin(1) = m.x(1); 
end