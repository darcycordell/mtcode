function [m]=load_model_wsinv(name)
%Function which loads a WSINV3D model file and puts it
%into a consistent format
%
% Usage: [m] = load_model_wsinv(name)
%
% "m" is the model structure
% "name" is the WSINV model (usually _model.XX or *.model) filename string
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
% -The assumption is that the WSINV model is referenced to the northwest top
% corner of the model. Is this always the case?

% -The calculation of the number of padding cells assumes a uniform mesh in
% the area of interest in the x and y directions. If some of the
% non-padding cells are different sizes, it will cause problems.
%
%
%
disp('Loading WSINV Model File')

m.name = name;
m.niter = ''; %Number of iterations left blank and can be added later

fid=fopen(name,'r');
fgetl(fid);
line=fscanf(fid,'%f %f %f %f',4);
m.nx=line(1);       %number of blocks in NS
m.ny=line(2);       %number of y blocks in EW
m.nz=line(3);     %number of z slices
index_rho_values=line(4); % resistivity index
    %If index_rho_values = 1, then the model files contains resistivity
    %indices (this is usually the inversion model input format from wsinv).
    %In other words, this is the initial model file.
    
    %If index_rho_values = 0, then the model file contains actual
    %resistivity values in a long column (this is the inversion model
    %output format from wsinv).
    
%Get cell thicknesses
m.dx=fscanf(fid,'%f',m.nx);  %block size
m.dy=fscanf(fid,'%f',m.ny);
m.dz=fscanf(fid,'%f',m.nz);

%Vector of cell edge locations starting from zero
tmp = [0; cumsum(m.dy)];
m.y = tmp-(max(tmp)/2);

tmp = [0; cumsum(m.dx)];
m.x = tmp-(max(tmp)/2);

m.origin = [m.x(1) m.y(1) 0]; %Currently, WSINV does not include elevations so the top of the model is set to 0.


tmp = cumsum(m.dz);
m.z = m.origin(3)+[0;tmp];

if index_rho_values % refers to resistivity indices
            % we are assuming this is the initial model (output uses absolute resistivities)
    
    fgetl(fid); % junk
    rho_ind = str2num(fgetl(fid)); % the resistivity values used
    m.A = zeros(m.nx,m.ny,m.nz);
    
    while 1
        layer = fscanf(fid,'%f',2); % read z start and end
        tmp = fscanf(fid,'%f',m.nx*m.ny); % read rho values for layer
        for ix = 1:m.nx*m.ny
            tmp(ix) = rho_ind(tmp(ix)); % replace index with actual rho values
        end
        rho_layer = reshape(tmp,m.nx,m.ny);    
        for il = layer(1):layer(end)
            m.A(:,:,il) = rho_layer;
        end
        if layer(end) == m.nz
            break
        end
    end
 
else % file uses actual resistivity values; assuming reading output model file

    m.A=reshape(fscanf(fid,'%f',m.nx*m.ny*m.nz),m.nx,m.ny,m.nz);      %z is the depth slices

end

fclose(fid);

%Elevations above topography surface are NaN
m.A(m.A>10^15) = NaN;

%WSINV is referenced from the north-west corner. To make it consistent with the southwest-origin
%format (common to many algorithms such as ModEM), it is necessary to
%flip the North-South dimension of the model.
m.A = m.A(end:-1:1,:,:);

%Find cell center midpoints
m.cx = (m.x(1:end-1)+m.x(2:end))/2;
m.cy = (m.y(1:end-1)+m.y(2:end))/2;
m.cz = (m.z(1:end-1)+m.z(2:end))/2;

%Build meshgrid of cell centers
[m.X,m.Y]=meshgrid(m.cy,m.cx);


m.nx = length(m.cx);
m.ny = length(m.cy);
m.nz = length(m.cz);

%Determine topography surface (should always be zeros since WSINV does not
%include topography)
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

%Determine padding
m.npad(1) = (m.nx - length(m.dx(m.dx==m.dx(round(m.nx/2)))))/2;
m.npad(2) = (m.ny - length(m.dy(m.dy==m.dy(round(m.ny/2)))))/2;


end % END function