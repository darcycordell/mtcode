function plot_slice_multi(m,d)
% Function which plots multiple NS, EW or horizontal model slices in 3-D view
%
% Usage: plot_slice_multi(m,d)
%
% Inputs: "m" is a standard model structure
%       "d" is a standard data structure (OPTIONAL)
%
% You can set the vertical exaggeration in user_defaults.m to separate the
% horizontal slices. See "u.ve" in user_defaults.m

u = user_defaults;
[L] = load_geoboundary_file_list;
close all

%If d structure does not exist then set the variable to all NaN
if ~exist('d','var')
    [d] = make_nan_data;
end

irun = 1;
set_figure_size(1);
while irun

main_menu = menu('Choose plane for slices','North-South','East-West','Horizontal layers','Return');
if main_menu == 4 || main_menu == 0
    irun = 0;
    break
end
prompt = 'Enter distances (or elevations b.s.l.) to slice, separated by a space';
dlg_title = 'Enter slice coordinates (km)';
def = {'-2 0 2'};
num_lines = 1;
inp = inputdlg(prompt,dlg_title,num_lines,def);
slice_num = str2num(char(inp));

% make 3-D arrays for X,Y,Z
X = repmat(m.X,[1,1,m.nz]);
Y = repmat(m.Y,[1,1,m.nz]);
tz = permute(m.z,[3 2 1]);
tz(:,:,end)=[];
Z=repmat(tz,[m.nx,m.ny,1]);

% trim padding cells
X = X(m.npad(1)+1:m.nx-m.npad(1),m.npad(2)+1:m.ny-m.npad(2),:)/1000;
Y = Y(m.npad(1)+1:m.nx-m.npad(1),m.npad(2)+1:m.ny-m.npad(2),:)/1000;
Z = Z(m.npad(1)+1:m.nx-m.npad(1),m.npad(2)+1:m.ny-m.npad(2),:)/1000;
M = m.A(m.npad(1)+1:m.nx-m.npad(1),m.npad(2)+1:m.ny-m.npad(2),:);

if main_menu == 1 % plotting North-South slices
    slice(X,Y,Z,log10(M),slice_num,[],[]);    
elseif main_menu == 2 % plotting East-West slices
    slice(X,Y,Z,log10(M),[],slice_num,[])
elseif main_menu == 3 % plotting horizontal layers
    slice(X,Y,Z,log10(M),[],[],slice_num)
end

hold on
% slice(X,Y,Z,log10(M),[],[],0) % testing


shading flat

colormap(u.cmap)
caxis(u.colim)
add_rho_colorbar(u);

set(gca,'dataaspectratio',[1 1 1/u.ve])
xlabel('Easting (km)')
ylabel('Northing (km)')
zlabel('Elevation b.s.l. (km)')

% axis([sort([m.cy(m.npad(2)+1) m.cy(m.ny-m.npad(2))]) sort([m.cx(m.npad(1)+1) m.cx(m.nx-m.npad(1))]) -u.zmax*1000 -u.zmin*1000]/1000)
y_lim = sort([m.cy(m.npad(2)+1) m.cy(m.ny-m.npad(2))])./1000;
x_lim = sort([m.cx(m.npad(1)+1) m.cx(m.nx-m.npad(1))])./1000;

if main_menu == 1 || main_menu == 2
    axis([y_lim x_lim u.zmin u.zmax])
else
    if length(slice_num)==1
        axis([y_lim x_lim u.zmin u.zmax])
    else
        axis([y_lim x_lim min(slice_num) max(slice_num)])
    end
end
set(gca,'zdir','reverse')

plot3(d.y/1000,d.x/1000,ones(size(d.z)).*min(slice_num),'k.','markerfacecolor','k') % project stations onto highest slice
% plot3(d.y/1000,d.x/1000,-(d.z-100)/1000,'kv','markerfacecolor','k') % switched x and y here instead of in the meshgrid
% plot3(d.x,d.y,d.z,'kv')

d.origin(3) = min(slice_num).*1000;
plot_geoboundaries(L,d.origin,d.origin(3));
% plot_geoboundaries(L,d.origin,d.z);


if main_menu == 1
    recz = [u.zmin u.zmax u.zmax u.zmin u.zmin];
    recx = [x_lim(1) x_lim(1) x_lim(2) x_lim(2) x_lim(1)];
    
    for i = 1:length(slice_num)
        recy = ones(size(recx)).*slice_num(i);
        plot3(recy,recx,recz,'k-')        
    end
    
elseif main_menu == 2
    recz = [u.zmin u.zmax u.zmax u.zmin u.zmin];
    recy = [y_lim(1) y_lim(1) y_lim(2) y_lim(2) y_lim(1)];
    
    for i = 1:length(slice_num)
        recx = ones(size(recy)).*slice_num(i);
        plot3(recy,recx,recz,'k-')        
    end   
    
else
    recy = [y_lim(1) y_lim(2) y_lim(2) y_lim(1) y_lim(1)];
    recx = [x_lim(2) x_lim(2) x_lim(1) x_lim(1) x_lim(2)];
    for i = 1:length(slice_num)
        recz = ones(size(recy)).*slice_num(i);
        plot3(recy,recx,recz,'k-')        
    end
end

camproj('perspective')
box on
set(gca,'BoxStyle','full')
set(gca,'fontsize',20)

set(gca,'view',[-29 35.5])


% figure % for rotating, still need to choose the best way to input the slice planes
% for k = 1:length(slice_num)
%     hs = slice(X,Y,Z,log10(M),slice_num(k),[],[]);
%     rotate(hs,[0 0 1],30)
%     xd = hs.XData; yd = hs.YData; zd = hs.ZData;
% %     delete(hs)
%     slice(X,Y,Z,log10(M),xd,yd,zd)
%     hold on
% end
end

end