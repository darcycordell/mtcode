function [p1] = bostick_penetration(fig_number,Z,Zvar,T,azimuth,bostick_grid,bostick_int_z,bostick_depth,stn_name,lat,long,elev,maplim,ns,profile_dir,geo_bound_file)
% 
% Plots a map of the Bostick resistivity at a set of depths for a grid
% of MT stations.  This is a simple version, calculating the Bostick
% resistivity from the Berdichevsky aveage.

% bostick_gridding2.m is more complex and calculates the maximum Bostick
% resistivity, and the azimuth in which that occurs. The rho_max direction
% is plotted on the map.

rad = 180./pi;      p1 = 1;   prof_azim = mean(azimuth);

mkdir bostick_penetration

% Process co-ords
x = long; y=lat;   
% x_mean = mean(x);  y_mean = mean(y); 
% x = cos(y_mean/rad)*111*(x-x_mean);    y = 111*(y-y_mean); % km
c = cosd(prof_azim);      s = sind(prof_azim);     R = [ c, -s ; s, c];

for is=1:ns
    loc = [x(is),y(is)];
    loc = R*loc';
    x_rot(is) = loc(1);y_rot(is)=loc(2);
end

% put stations in order on monotonically increasing x_rot
[x_sort,index] = sort(x_rot,'ascend');

   
%==========================================================================
% calculate the Bostick resistivity from the Berdichevsky average
for is =1:ns    
    i_sort = index(is);  
    y_sort(is)=y_rot(i_sort);
    Z1 = Z(:,:,:,i_sort);    Z1var = Zvar(:,:,:,i_sort);
    [rho_bostick_av(:,is),d_av(:,is),~,~] = bostick1(Z1,T);
end

% Plot map to check co-ords. If all dots on same position then no problem  
% x,y is the original one, x_rot, y_rot is rotated one, x_sort, y_sort is 
% sortted in order on monotonically increasing x_rot.
figure(100)
m_grid('tickdir','out','yaxislocation','left','xlabeldir','end','ticklen',.02); 
hold on; 
m_plot(x,y,'o');
m_plot(x_rot,y_rot,'+');
m_plot(x_sort,y_sort,'r.');
title('COORDS CHECK - press any key to continue');

pause (2);
close(100);

%==========================================================================
% Set up the grid  (degree) for lateral interpolation (griddata)
% x_gridding=x_sort;     y_gridding=y_sort;
% xg = (min(x_gridding):bostick_grid(1):max(x_gridding));               
% yg = (min(y_gridding):bostick_grid(2):max(y_gridding));
xg = (maplim(1):bostick_grid(1):maplim(2));               
yg = (maplim(3):bostick_grid(2):maplim(4));
[X,Y] = meshgrid(xg,yg); 

%==========================================================================
for i=1:size(d_av,1)
    d_gridding = log10(d_av(i,:)/1000);
    x_gridding=x_sort;     y_gridding=y_sort;
    % delete the data corresponding to NAN in bostick_pseudo
    x_gridding(isnan(d_gridding))=[];
    y_gridding(isnan(d_gridding))=[];
    d_gridding(isnan(d_gridding))=[];
    if length(x_gridding)<=5
        d_grid_av(:,:,i) = nan(size(X));
    else
%         bostick_grid_av(:,:,i) = griddata(x_gridding,y_gridding,bostick_rho_gridding,X,Y,'nearest');
%         F = TriScatteredInterp(x_gridding',y_gridding',d_gridding','natural');
        F = scatteredInterpolant(x_gridding',y_gridding',d_gridding','natural','linear');
        % You can use natural or nearest as the interpolation method, and
        % nearest is closer to the real data but rougher than the natural
        % method.
        d_grid_av(:,:,i) = F (X,Y);
    end
end
%==========================================================================
% Plot the penetration depth of Bostick results at different periods

for i=1:length(T) 
    figure(fig_number);
    m_grid('tickdir','out','yaxislocation','left','xlabeldir','end','ticklen',.02); 
    hold on;
    m_pcolor(xg,yg,d_grid_av(:,:,i));  
    m_plot(x,y,'k.');
    plot_geo_bound_file_m_map(geo_bound_file);
    hold off;
    
    colormap(flipud(jet));    colorbar vertical;
    shading flat;    
%     axis([min(xg),max(xg),min(yg),max(yg)]);
    caxis([0,3]);
%     caxis auto;
    title(['Depth Penetration of Average Resistivity Bostick T = ',...
        num2str(T(i),'%6.4f'),'  second'])
%     xlabel('longitude'); ylabel('latitude'); 
    hold on; 
    

    % output the figure
    eval(['print -djpeg100 bostick_penetration\bostick_penetration_depth_',...
        num2str(T(i),'%6.4f'),'_second.jpg']);
    
    pause (2)
    close(fig_number)

end 
end