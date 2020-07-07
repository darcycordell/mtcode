function [m] = load_model_2DNLCG(m_file,sitefile)
%Function which loads Randy Mackie's 2D NLCG Inversion format model files
%and puts it into a consistent format in a "m" (model) structure.
%
% The model structure is standardized and is compatible with both 2D and 3D
% models. For the 2D case, the x-direction is said to be perpendicular to
% the profile and the y-direction is along profile. As such, the
% x-direction is only 1 cell thick, 1 m wide.
%
%The model structure includes:
% A = Model resistivity values
% dx = Model cell widths in the off-profile direction (Set to "1" for 2D)
% dy = Model cell widths in the east-west direction in meters
% dz = Model cell widths in depth in meters
% x = South cell edge vertex locations referenced to center of mesh
% y = West cell edge vertex locations referenced to center of mesh
% z = Top cell edge vertex locations in depth below sea level
% X = Meshgrid of cell edges in north-south direction
% Y = Meshgrid of cell edges in east-west direction
% origin = The top, southwest corner of the model
    %Note: The elevation of the top of the model must be referend to meters below
    %sea level (i.e. elevations which are above sea level should be
    %negative; if you model has topography then origin(3) < 0)
% cx = the cell centers in the north-south direction. Since the resistivity
% values are defined at cell centers, it is better to plot the resistivity
% values aligned with cell centers
% cy = cell centers in the east-west direction.
% X = Meshgrid of east-west edges (north-south map view)
% Y = Meshgrid of north-south edges (east-west map view)
%
% Note: Because "dy" are cell thickness and "y" are cell edges, the length
% of dy is N and the length of y is (N+1).
%
% To be added:----------------------------------------------------------
%
% -How to get the origin of the mesh? Relative to what in lat-long space?
% -Currently the number of padding cells is set to zero, not sure how to
% count where padding cells being and end on a WinGLink mesh.


fid=fopen(m_file,'r');

disp('Loading 2D NLCG Model File')

[~,~,ref] = load_sitefile(sitefile,1);

tmp=fscanf(fid,'%i',2);
m.nx = 1; %For 2D, the number of cells in the x-direction is 1.
m.ny=tmp(1);
m.nz=tmp(2);

m.origin = [0,0,0]; %how to get origin of mesh? Relative to what?

m.dx = 1; %The thickness in the x-direction (perpendicular to profile) is set to 1 m
m.dy=fscanf(fid,'%f',m.ny);
m.dz=fscanf(fid,'%f',m.nz);

tmp = fscanf(fid,'%i',1);

tmp=fscanf(fid,'%f',m.ny*m.nz*m.nx);

tline = fgetl(fid);
tline = fgetl(fid);
while 1    
    if ~isempty(strfind(tline,'#'))
        tline = fgetl(fid);
    else
        break
    end
end
m.z_ind = str2num(tline); % get the ksurf values at the bottom of the model file

fclose(fid);

m.A=reshape(tmp,m.nx,m.ny,m.nz);


%Vector of cell edge locations starting from zero
tmp = cumsum(m.dy);
tmp = [0;tmp];
m.y = tmp+ref;

tmp = cumsum(m.dx);
tmp = [0;tmp];
m.x = tmp-(max(tmp)/2);

tmp = cumsum(m.dz);
m.z = m.origin(3)+[0;tmp]; %how to get true elevation?

m.A(m.A>=10^10) = NaN;

%Find cell center midpoints
m.cx = (m.x(1:end-1)+m.x(2:end))/2;
m.cy = (m.y(1:end-1)+m.y(2:end))/2;
m.cz = (m.z(1:end-1)+m.z(2:end))/2;

% [m.Y,m.X]=meshgrid(m.cy,m.cx); 
[m.Y,m.X,m.Z]=meshgrid(m.cy,m.cx,m.cz); % some indices mixed up? maybe fixed 9/20/2019
% m.Z = zeros(m.nx,m.ny);
% for i = 1:m.ny
%     for j = 1:m.nz
%         
%         ind = find(isnan(squeeze(m.A(i,j,:))),1,'last');
%         
%         if isempty(ind)
%             ind = 0;
%         end
%         
%         m.Z(i,j) = m.cz(ind+1);
%         j
%     end
%     
% end

m.npad = [0, 0]; %How to define "padding cells" for a WinGLink mesh?




end