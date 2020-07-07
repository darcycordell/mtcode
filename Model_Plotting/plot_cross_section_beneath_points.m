function [x,y,z,rho]=plot_cross_section_beneath_points(m,d)
% Function which plots a cross-section fence diagram between selected MT sites
% Points can be selected by clicking on a map or by supplying a text file
%
% Usage: [x,y,z,rho] = plot_cross_section_beneath_points(m,d)
%
% Inputs:
% "m" is the model structure, "d" is the data structure
%
% Outputs:
% x,y,z are the NS, EW and depth locations of the section
% rho is the matrix of resistivity values associated with the section
%
% If a text file is supplied, it must be either:
% two columns: (1) latitudes and (2) longitudes,
% a single column list of station indices
%
u = user_defaults;
menu_1 = menu('Get points from:','Text File','Select from Map');
close all
if menu_1 == 1
    curdir = pwd;
    [profile_points_file, filepath]=uigetfile({'*.txt'},'Choose text file which contains longitude and latitude points OR station indices'); if profile_points_file == 0; return; end
    cd(filepath)

    % Define profile of chosen stations and place nseg points between each pair 
    profile_points=load(profile_points_file);
    fclose all;
    cd(curdir)
    
    if size(profile_points,2)==1
        
        x = d.x(profile_points)/1000;
        y = d.y(profile_points)/1000;
        
    elseif size(profile_points,2)==2
    
        [y,x] = geo2utm(profile_points(:,2),profile_points(:,1),d.origin(2),d.origin(1));
        y = (y-500000)/1000;
        x = x/1000;
        
    else
        disp('***Error: Something is wrong with your text file specifying locations to extract profile along')
        disp('See S3 documentation. The text file should either contain one column of station indices OR')
        disp('Two columns: latitude followed by longitude of points from which to extract profile')
        return
    end
    
    plot_slice(m,round(m.nz/2),d)
    
    ind = zeros(length(x),1);
    for i = 1:length(x)
        dist = sqrt((x(i)-d.x/1000).^2+(y(i)-d.y/1000).^2);
        [~,ind(i)] = min(dist);

        x(i) = d.x(ind(i))/1000;
        y(i) = d.y(ind(i))/1000;
        
    end

elseif menu_1 == 2
    
    plot_slice(m,round(m.nz/2),d)
    i=0; x=0; y=0; ind=0;
    while 1
        [yp,xp,click] = ginput(1); 
        if click~=3
            i=i+1;
            dist = sqrt((xp-d.x/1000).^2+(yp-d.y/1000).^2);
            [~,ind(i)] = min(dist);

            x(i) = d.x(ind(i))/1000;
            y(i) = d.y(ind(i))/1000;
            
            plot(y(i),x(i),'.r','MarkerSize',14)

        else
            break
        end
        
        
    end
    
end

plot(y,x,'.-r','MarkerSize',14,'LineWidth',3)

%%%%%%%%%%%%%%%%%%
nseg=10; %  set to 1 for just the values at each station
%%%%%%%%%%%%%%%%%%

ipts=1; % Initialize counter
np = length(x);
px = zeros((length(x)-1)*nseg,1);
py = zeros((length(x)-1)*nseg,1);
for isp =1:np-1   % Start loop over list of points
  dx = x(isp+1)-x(isp);   
  dy = y(isp+1)-y(isp);
  for i=1:nseg
    px(ipts)= x(isp) + (i-1)*dx/nseg; 
    py(ipts)= y(isp) + (i-1)*dy/nseg;
    ipts=ipts+1;
  end
end

for i=1:length(px)
    plot(py(i),px(i),'r*');hold on
end

dist = zeros(length(px),1);
for i = 2:length(px)
    dist(i) = sqrt((px(i)-px(i-1)).^2+(py(i)-py(i-1)).^2);
end

dist = cumsum(dist);
    

zind = 1:nearestpoint(u.zmax*1000,m.cz);
% Now interpolate, layer by layer
res_plot = zeros(length(px),length(zind));
for iz=zind
 slice = m.A(:,:,iz);
 res_plot(:,iz) =interp2(m.X/1000,m.Y/1000,slice,py,px);
end


figure(2)  
pcolor(dist,m.cz(zind)/1000,log10(res_plot)'); axis ij; shading flat; hold on
colormap(u.cmap); caxis(u.colim);
add_rho_colorbar(u)

xlabel('Distance along Profile (km)');
ylabel('Depth (km)');
axis([min(dist) max(dist) u.zmin u.zmax]);
set(gca,'Layer','top')

tstr = [d.site{ind(1)},'_to_',d.site{ind(end)},'_',num2str(numel(ind)),'_stations'];
title(strrep(['Slice Through Model at Sites ',tstr],'_','\_'));
set(gca,'DataAspectRatio',[u.ve 1 1])

plot([dist(1:nseg:end); dist(end)],d.z(ind)/1000 -u.zoff,'vk','MarkerFaceColor','k')

print_figure(['diagonal_slices_',m.niter],['slice_through_sites_',tstr]) 

figure(1)
print_figure(['diagonal_slices_',m.niter],['slice_through_sites',tstr,'_MAP']) 

x = px;
y = py;
z = m.cz(zind)/1000;
rho = res_plot(:,zind);



end
