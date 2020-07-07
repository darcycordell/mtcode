function plot_model_3D(m,d)
%
% Function which renders MT models in 3D as cell faces and vertices.
%
% Usage: plot_model_3D(m,d)
%
% Inputs: "m" is a standard model structure
%       "d" is a standard data structure
%
%

u = user_defaults;
[L] = load_geoboundary_file_list;
close all

z_lim = [u.zmin u.zmax];

plot_stations = 1; % 1 to plot stations, 0 to skip

% should not need to edit below here
f1 = 1:4; % vertices to connect for 1 face of a rectangular cell
f5 = [1 2 3 4; 1 2 6 5; 2 3 7 6; 3 4 8 7; 1 4 8 5]; % vertices to connect for 5 faces (except bottom) of a rectangular cell

[mx,my,mz] = meshgrid(m.x./1000,m.y./1000,(m.z)./1000);
mx = permute(mx,[2 1 3]);
my = permute(my,[2 1 3]);
mz = permute(mz,[2 1 3]);
A = log10(m.A);

set_figure_size(1);

while 1
lim_menu = menu('3D Plotting Options','Plot with default x and y limits','Enter custom x and y limits','Return');

if lim_menu == 1 % find default limits
    xind = m.npad(1)+1:m.nx-m.npad(1);
    yind = m.npad(2)+1:m.ny-m.npad(2);
    xn_ind = xind(end);
    xs_ind = xind(1);   
    ye_ind = yind(end);
    yw_ind = yind(1);   
elseif lim_menu == 2 % enter custom x, y limits and use z limits from user_defaults  
    prompt = {'x min (km)';'x max (km)';'y min (km)';'y max (km)'};
    dlg_title = 'Plot limits';
    def = {'-20'; '20'; '-20'; '20'};
    num_lines = 1;
    inp = inputdlg(prompt,dlg_title,num_lines,def);
    
    xn_ind = nearestpoint(str2double(inp{2}),unique(mx(:,:,1))); % find min and max indices to plot in each dim
    xs_ind = nearestpoint(str2double(inp{1}),unique(mx(:,:,1)));
    ye_ind = nearestpoint(str2double(inp{4}),unique(my(:,:,1)));
    yw_ind = nearestpoint(str2double(inp{3}),unique(my(:,:,1)));
else
    close all
    return
end
clf

zb_ind = nearestpoint(u.zmax,unique(mz(1,1,:)));
zt_ind = nearestpoint(u.zmin,unique(mz(1,1,:)));
%% top face and air layers
topobot_ind = zeros(xn_ind-xs_ind+1,ye_ind-yw_ind+1);
for iz = 1:zb_ind % find lowest points in topo
    tmp = isnan(A(xs_ind:xn_ind,yw_ind:ye_ind,iz));
    topobot_ind = topobot_ind + tmp;
    if sum(sum(tmp)) ==0 % if entirely below ground
        break
    end
end
% max(max(topobot_ind)) + 1 is the first layer with no air cells in chosen plot area
% min(min(topobot_ind)) + 1 is the first layer with an earth cell

% 1. Check layer to plot. 
% a) If in air, loop from layer to first layer with earth cell. 
% In each step find the cells to be plotted and combine into one patch object.
% b) If below ground, just plot that one slice.

if zt_ind > max(max(topobot_ind)) % if layer to be plotted has no air cells
    tic
    
    verts = zeros((ye_ind-yw_ind)*(xn_ind-xs_ind)*4,3); % 4 vertices per cell
    faces = zeros((ye_ind-yw_ind)*(xn_ind-xs_ind),4); % just one face per cell
    colors = zeros((ye_ind-yw_ind)*(xn_ind-xs_ind),1);
    
    count = 1;
    for iew = yw_ind:ye_ind-1
        for ins = xs_ind:xn_ind-1
        v = [[mx(ins,iew,zt_ind) my(ins,iew,zt_ind) mz(ins,iew,zt_ind)];...
        [mx(ins+1,iew,zt_ind) my(ins+1,iew,zt_ind) mz(ins+1,iew,zt_ind)];...
        [mx(ins+1,iew+1,zt_ind) my(ins+1,iew+1,zt_ind) mz(ins+1,iew+1,zt_ind)];...
        [mx(ins,iew+1,zt_ind) my(ins,iew+1,zt_ind) mz(ins,iew+1,zt_ind)]];
               
        verts((count*4)-3:count*4,:) = v; % done for count = 1
        faces(count,:) = f1 + (count-1)*4;
        colors(count,:) =  A(ins,iew,zt_ind); 
        
        count = count + 1;
        end
    end    
    
    disp(['Layer ',num2str(zt_ind)])
    toc

    patch('vertices',[verts(:,2) verts(:,1) verts(:,3)] ,'faces',faces,'cdata',colors,'facecolor','flat') % x and y are swapped here
    hold on
    
else % if layer to be plotted contains air cells
           
    for iz = zt_ind:(max(max(topobot_ind)) + 1)
%         tic
        chk = sum(sum(~isnan(A(xs_ind:xn_ind,yw_ind:ye_ind,iz)))); % if 0, then this layer is all air cells and no plotting needed
        if chk > 0
            verts = zeros((ye_ind-yw_ind)*(xn_ind-xs_ind)*8,3); % going to be 8 vertices per cell, each vertex is a xyz trio
            faces = zeros((ye_ind-yw_ind)*(xn_ind-xs_ind)*5,4); % going to be 5 faces per cell, 4 vertices per face
            colors = zeros((ye_ind-yw_ind)*(xn_ind-xs_ind)*5,1); % one color per face
            count = 1;
            for iew = yw_ind:ye_ind-1
                for ins = xs_ind:xn_ind-1
                    v = [[mx(ins,iew,iz) my(ins,iew,iz) mz(ins,iew,iz)];... % these are the 8 vertices for each cell
                    [mx(ins+1,iew,iz) my(ins+1,iew,iz) mz(ins+1,iew,iz)];...
                    [mx(ins+1,iew+1,iz) my(ins+1,iew+1,iz) mz(ins+1,iew+1,iz)];...
                    [mx(ins,iew+1,iz) my(ins,iew+1,iz) mz(ins,iew+1,iz)]                     
                    [mx(ins,iew,iz+1) my(ins,iew,iz+1) mz(ins,iew,iz+1)];...
                    [mx(ins+1,iew,iz+1) my(ins+1,iew,iz+1) mz(ins+1,iew,iz+1)];...
                    [mx(ins+1,iew+1,iz+1) my(ins+1,iew+1,iz+1) mz(ins+1,iew+1,iz+1)];...
                    [mx(ins,iew+1,iz+1) my(ins,iew+1,iz+1) mz(ins,iew+1,iz+1)]];

                    verts( (8*count)-7:8*count ,:) = v; 
                    faces(((count-1)*5)+1:count*5,:) = f5+((count-1)*8); 
                    colors(((count-1)*5)+1:count*5,:) =  A(ins,iew,iz); 

                    count = count + 1;
                end
            end
            % find and remove vertices corresponding to air cells. these are in sets of 8
            tmp = find(isnan(colors(1:5:end,:))).*8; % since one color per face, need to find colors in blocks of 5, then multiply by 8 to include all verts 
            verts(([tmp; tmp-1; tmp-2; tmp-3; tmp-4; tmp-5; tmp-6; tmp-7]),:) = []; % need to remove all 8 verts per cell at once, thus the substractions
            faces = faces(1:size(verts,1)/8*5,:); % truncate to match 5 faces per 8 verts
            colors(isnan(colors)) = [];

            patch('vertices',[verts(:,2) verts(:,1) verts(:,3)] ,'faces',faces,'cdata',colors,'facecolor','flat') % x and y are swapped here
            hold on        

        end
%         disp(['Layer ',num2str(iz)])
%         toc
    end                          
    
end % end plotting top faces
%%
% north face
% tic

verts = zeros((ye_ind-yw_ind)*(zb_ind-zt_ind)*4,3); % 4 vertices per cell
faces = zeros((ye_ind-yw_ind)*(zb_ind-zt_ind),4); % just one face per cell
colors = zeros((ye_ind-yw_ind)*(zb_ind-zt_ind),1);
count = 1;

for ie = yw_ind:ye_ind-1
    for iz = zt_ind:zb_ind-1
        v = [[mx(xn_ind,ie,iz) my(xn_ind,ie,iz) mz(xn_ind,ie,iz)];...
        [mx(xn_ind,ie+1,iz) my(xn_ind,ie+1,iz) mz(xn_ind,ie+1,iz)];...
        [mx(xn_ind,ie+1,iz+1) my(xn_ind,ie+1,iz+1) mz(xn_ind,ie+1,iz+1)];...
        [mx(xn_ind,ie,iz+1) my(xn_ind,ie,iz+1) mz(xn_ind,ie,iz+1)]];
    
        verts((count*4)-3:count*4,:) = v;
        faces(count,:) = f1 + (count-1)*4;
        colors(count,:) =  A(xn_ind-1,ie,iz); % because of indexing the SW corner of each cell, we need to plot the north face of each cell by indexing ix-1
        
        count = count + 1;                              
    end
end

verts(([find(isnan(colors)).*4; (find(isnan(colors)).*4)-1; (find(isnan(colors)).*4)-2; (find(isnan(colors)).*4)-3 ]),:) = [];
faces = faces(1:length(verts)/4,:); % truncate to match number of vertices / 4 
colors(isnan(colors)) = [];

% disp('North Face (c)')
% toc

patch('vertices',[verts(:,2) verts(:,1) verts(:,3)] ,'faces',faces,'cdata',colors,'facecolor','flat') % x and y are swapped here
hold on

% south face
% tic

verts = zeros((ye_ind-yw_ind)*(zb_ind-zt_ind)*4,3); % 4 vertices per cell
faces = zeros((ye_ind-yw_ind)*(zb_ind-zt_ind),4); % just one face per cell
colors = zeros((ye_ind-yw_ind)*(zb_ind-zt_ind),1);
count = 1;

for iw = yw_ind:ye_ind-1
    for iz = zt_ind:zb_ind-1
        v = [[mx(xs_ind,iw,iz) my(xs_ind,iw,iz) mz(xs_ind,iw,iz)];...
        [mx(xs_ind,iw+1,iz) my(xs_ind,iw+1,iz) mz(xs_ind,iw+1,iz)];...
        [mx(xs_ind,iw+1,iz+1) my(xs_ind,iw+1,iz+1) mz(xs_ind,iw+1,iz+1)];...
        [mx(xs_ind,iw,iz+1) my(xs_ind,iw,iz+1) mz(xs_ind,iw,iz+1)]];
               
        verts((count*4)-3:count*4,:) = v;
        faces(count,:) = f1 + (count-1)*4;
        colors(count,:) =  A(xs_ind,iw,iz); 
        
        count = count + 1;                
    end
end

verts(([find(isnan(colors)).*4; (find(isnan(colors)).*4)-1; (find(isnan(colors)).*4)-2; (find(isnan(colors)).*4)-3 ]),:) = [];
faces = faces(1:length(verts)/4,:); % truncate to match number of vertices / 4 
colors(isnan(colors)) = [];

% disp('South Face')
% toc

patch('vertices',[verts(:,2) verts(:,1) verts(:,3)] ,'faces',faces,'cdata',colors,'facecolor','flat') % x and y are swapped here

% west face
% tic

verts = zeros((xs_ind-xn_ind)*(zb_ind-zt_ind)*4,3); % 4 vertices per cell
faces = zeros((xs_ind-xn_ind)*(zb_ind-zt_ind),4); % just one face per cell
colors = zeros((xs_ind-xn_ind)*(zb_ind-zt_ind),1);
count = 1;

for in = xs_ind:xn_ind-1
    for iz = zt_ind:zb_ind-1
        v = [[mx(in,yw_ind,iz) my(in,yw_ind,iz) mz(in,yw_ind,iz)];...
        [mx(in+1,yw_ind,iz) my(in+1,yw_ind,iz) mz(in+1,yw_ind,iz)];...
        [mx(in+1,yw_ind,iz+1) my(in+1,yw_ind,iz+1) mz(in+1,yw_ind,iz+1)];...
        [mx(in,yw_ind,iz+1) my(in,yw_ind,iz+1) mz(in,yw_ind,iz+1)]];
                  
        verts((count*4)-3:count*4,:) = v;
        faces(count,:) = f1 + (count-1)*4;
        colors(count,:) =  A(in,yw_ind,iz); 
        
        count = count + 1;            
    end
end

verts(([find(isnan(colors)).*4; (find(isnan(colors)).*4)-1; (find(isnan(colors)).*4)-2; (find(isnan(colors)).*4)-3 ]),:) = [];
faces = faces(1:length(verts)/4,:); % truncate to match number of vertices / 4 
colors(isnan(colors)) = [];

% disp('West Face')
% toc

patch('vertices',[verts(:,2) verts(:,1) verts(:,3)] ,'faces',faces,'cdata',colors,'facecolor','flat') % x and y are swapped here

% east face
% tic

verts = zeros((xs_ind-xn_ind)*(zb_ind-zt_ind)*4,3); % 4 vertices per cell
faces = zeros((xs_ind-xn_ind)*(zb_ind-zt_ind),4); % just one face per cell
colors = zeros((xs_ind-xn_ind)*(zb_ind-zt_ind),1);
count = 1;

for is = xs_ind:xn_ind-1
    for iz = zt_ind:zb_ind-1
        v = [[mx(is,ye_ind,iz) my(is,ye_ind,iz) mz(is,ye_ind,iz)];...
        [mx(is+1,ye_ind,iz) my(is+1,ye_ind,iz) mz(is+1,ye_ind,iz)];...
        [mx(is+1,ye_ind,iz+1) my(is+1,ye_ind,iz+1) mz(is+1,ye_ind,iz+1)];...
        [mx(is,ye_ind,iz+1) my(is,ye_ind,iz+1) mz(is,ye_ind,iz+1)]];    
        
        verts((count*4)-3:count*4,:) = v;
        faces(count,:) = f1 + (count-1)*4;
        colors(count,:) =  A(is,ye_ind-1,iz);  % iy-1 because we index the SW corner of each cell
        
        count = count + 1;   
    end
end

verts(([find(isnan(colors)).*4; (find(isnan(colors)).*4)-1; (find(isnan(colors)).*4)-2; (find(isnan(colors)).*4)-3 ]),:) = [];
faces = faces(1:length(verts)/4,:); % truncate to match number of vertices / 4 
colors(isnan(colors)) = [];

% disp('East Face')
% toc

patch('vertices',[verts(:,2) verts(:,1) verts(:,3)] ,'faces',faces,'cdata',colors,'facecolor','flat') % x and y are swapped here

% bottom face 
% tic

verts = zeros((ye_ind-yw_ind)*(xn_ind-xs_ind)*4,3); % 4 vertices per cell
faces = zeros((ye_ind-yw_ind)*(xn_ind-xs_ind),4); % just one face per cell
colors = zeros((ye_ind-yw_ind)*(xn_ind-xs_ind),1);
count = 1;
for iew = yw_ind:ye_ind-1
    for ins = xs_ind:xn_ind-1
    v = [[mx(ins,iew,zb_ind) my(ins,iew,zb_ind) mz(ins,iew,zb_ind)];...
    [mx(ins+1,iew,zb_ind) my(ins+1,iew,zb_ind) mz(ins+1,iew,zb_ind)];...
    [mx(ins+1,iew+1,zb_ind) my(ins+1,iew+1,zb_ind) mz(ins+1,iew+1,zb_ind)];...
    [mx(ins,iew+1,zb_ind) my(ins,iew+1,zb_ind) mz(ins,iew+1,zb_ind)]];

    verts((count*4)-3:count*4,:) = v; % done for count = 1
    faces(count,:) = f1 + (count-1)*4;
    colors(count,:) =  A(ins,iew,zb_ind); 

    count = count + 1;
    end
end    

% disp('Bottom Face')
% toc

patch('vertices',[verts(:,2) verts(:,1) verts(:,3)] ,'faces',faces,'cdata',colors,'facecolor','flat') % x and y are swapped here

%% make the plot look nice
view(3)
%set(gca,'zdir','reverse')
% zlim(z_lim)

colormap(u.cmap); caxis(u.colim); 
add_rho_colorbar(u);

set(gca,'dataaspectratio',[1 1 1/u.ve])
set(gca,'zdir','reverse')
xlabel('Easting (km)')
ylabel('Northing (km)')
zlabel('Elevation b.s.l (km)')
zlabel('Elevation b.s.l. (km)')
axis([m.cy(yw_ind-1) m.cy(ye_ind+1) m.cx(xs_ind-1) m.cx(xn_ind+1) u.zmin*1000+500 u.zmax*1000-500]/1000)
box on


if plot_stations
    if ~unique(d.z) % if all elevations are zero
        scatter3(d.y/1000,d.x/1000,-d.z/1000,'kv','markerfacecolor','k')
    else
        znew = zeros(size(d.z));
        for i = 1:length(d.z)
            dif = (d.z(i)-m.z)./1000;
            [~,ind] = min(dif(dif > 0));
            znew(i) = m.z(ind);
        end
        scatter3(d.y/1000,d.x/1000,znew/1000,'kv','markerfacecolor','k')
    end
end

plot_geoboundaries(L,d.origin,d.z)

end % end while

end % end function