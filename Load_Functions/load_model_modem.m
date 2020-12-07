function [m] = load_model_modem(name)
%Function which loads a ModEM model file and puts it
%into a consistent format
%
% Usage: [m] = load_model_modem(name)
%
% "m" is the model structure
% "name" is the ModEM model (usually *.rho or *.model) filename string
%
% Embedded function: [x,y,z,rho,nzAir,type,origin,rotation] = read_mackie3d_model(fname,block)
%       This function is from Anna Kelbert
%
%The model structure includes:
% A = Model resistivity values
% dx = Model cell widths in the north-south direction in meters
% dy = Model cell widths in the east-west direction in meters
% dz = Model cell widths in depth in meters
% nx = number of cells in x (N-S) direction
% ny = number of cells in y (E-W) direction
% nz = number of cells in z direction
% x = South cell edge vertex locations referenced to center of mesh
% y = West cell edge vertex locations referenced to center of mesh
% z = Top cell edge vertex locations in depth below sea level
% origin = The top, southwest corner of the model
%    Note: The elevation of the top of the model should be referend to meters below
%    sea level (i.e. elevations which are above sea level should be
%    negative; if you model has topography then origin(3) < 0)
% cx = the cell centers in the north-south direction. Since the resistivity
%   values are defined at cell centers, it is better to plot the resistivity
%   values aligned with cell centers.
% cy = cell centers in the east-west direction.
% cz = cell centers in depth
% X = Meshgrid of cell centers in north-south direction
% Y = Meshgrid of cell centers in east-west direction
% Z = elevation topography surface (nx x ny)
% npad = number of padding cells
%
% Note: Because "dx" are cell thickness and "x" are cell edges, the length
% of dx and cx is N and the length of x is (N+1).
%
% To be added:----------------------------------------------------------
%
% Mesh rotations are not currently implemented
%
%
% Potential Bugs: ------------------------------------------------------
%
% -The calculation of the number of padding cells assumes a uniform mesh in
% the area of interest in the x and y directions. If some of the
% non-padding cells are different sizes, it will cause problems.
%
%
%------------------------------------------------------------------------

disp('Loading ModEM Model File')
%Read model file using old function
[m.dx,m.dy,m.dz,m.A,~,type,m.origin,~] = read_mackie3d_model(name,'true');

m.name = name;
m.niter = ''; %Number of iterations left blank and can be added later

%If model is logarithmic then adjust for that
if strcmp(type,'LOGE')
    m.A = exp(m.A);
end

%Vector of cell edge locations starting from zero
% tmp = [0; cumsum(m.dy)]; %Old way which does not incorporate origin
% m.y = tmp-(max(tmp)/2);
m.y = [0; cumsum(m.dy)] + m.origin(2);

%tmp = [0; cumsum(m.dx)]; %Old way which does not incorporate origin
%m.x = tmp-(max(tmp)/2);
m.x = [0; cumsum(m.dx)] + m.origin(1);

m.z = [0; cumsum(m.dz)] + m.origin(3);

%Elevations above topography surface are set to NaN
m.A(m.A>10^15) = NaN;

%Find cell center midpoints
m.cx = (m.x(1:end-1)+m.x(2:end))/2;
m.cy = (m.y(1:end-1)+m.y(2:end))/2;
m.cz = (m.z(1:end-1)+m.z(2:end))/2;

%Create a meshgrid of cell edges
[m.X,m.Y]=meshgrid(m.y,m.x);
[m.Xc,m.Yc] = meshgrid(m.cy,m.cx);


m.nx = length(m.cx);
m.ny = length(m.cy);
m.nz = length(m.cz);

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

%Determine padding
m.npad(1) = (m.nx - length(m.dx(m.dx==m.dx(round(m.nx/2)))))/2;
m.npad(2) = (m.ny - length(m.dy(m.dy==m.dy(round(m.ny/2)))))/2;



end %END MAIN

function [x,y,z,rho,nzAir,type,origin,rotation] = read_mackie3d_model(fname,block)
% reads a 3D resistivity model in Randie Mackie's format
% if 'block' is specified and true, no integer indices are read
% and the model is transposed (provides support for WS's models)
if nargin < 2
    block = 0;
end
fid = fopen(fname);
line = fgetl(fid);
while line(1)=='#' || line(2)=='#'
    line = fgetl(fid);
end
[n] = sscanf(line,'%d',[4 1]);
if findstr(line,'LOGE')
    type = 'LOGE';
else
    type = 'LINEAR';
end
nx =   n(1);   ny  =   n(2);   nz  =   n(3);
nzAir = n(4);
x   =   fscanf(fid,'%f',nx);
y   =   fscanf(fid,'%f',ny);
z   =   fscanf(fid,'%f',nz);
k = 1;
rho(1:nx,1:ny,1:nz) = 0;
if block
    for i = 1:nz
        tmp = fscanf(fid,'%f',[nx ny]);
        rho(:,:,i) = flipud(tmp);
    end
else
while max(k) < nz
    tmp =   fscanf(fid,'%s\n',1);
    k   =   sscanf(tmp,'%d');
    if length(k)<2; k = [k k]; end
    tmp =   fscanf(fid,'%f',[nx ny]);
    for i = k(1):k(2)
       rho(:,:,i) = tmp;
    end
end
end
origin = [0 0 0];
rotation = 0;
%DC Edits: fix so that read model actually outputs the origin and rotation
%values
while 1
line = fgetl(fid); 
if ~ischar(line); fclose(fid); return; end

nline = str2num(line);

if isempty(nline)==0
    if length(nline)==3
        origin = nline;
    elseif length(nline)==1
        rotation = nline;
    end
end
end
    
fclose(fid);

end
