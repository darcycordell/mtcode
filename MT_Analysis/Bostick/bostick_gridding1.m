function [p1] = bostick_gridding1(fig_number,Z,Zvar,T,azimuth,rlim,bostick_grid,bostick_int_z,bostick_depth,stn_name,lat,long,elev,maplim,ns,profile_dir,geo_bound_file)
% 
% Plots a map of the Bostick resistivity at a set of depths for a grid
% of MT stations.  This is a simple version, calculating the Bostick
% resistivity from the Berdichevsky aveage.

% bostick_gridding2.m is more complex and calculates the maximum Bostick
% resistivity, and the azimuth in which that occurs. The rho_max direction
% is plotted on the map.

rad = 180./pi;      p1 = 1;   prof_azim = mean(azimuth);

mkdir bostick

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
% Using interpolation to get the average bostick resistivity model
logdepth=log10(d_av);   % Array of all depths (all stations, frequencies)
d_interp=min(min(logdepth)):bostick_int_z:max(max(logdepth));   % Define depths onto which we interpolate
  
for is=1:ns
    temp1 = logdepth(:,is);          temp1(isnan(temp1))=[]; % Remove the NaN values from the vector
    temp2 = rho_bostick_av(:,is);    temp2(isnan(temp2))=[]; % Remove the NaN values from the vector
    if length(temp1) <=5
        bostick_pseudo_av(:,is)=nan(size(d_interp));
    else
        bostick_pseudo_av(:,is)=abs(interp1(temp1,temp2,d_interp,'cubic',NaN));
    end
end

clear temp1 temp2;
%==========================================================================
% Set up the grid  (degree) for lateral interpolation
% x_gridding=x_sort;     y_gridding=y_sort;
% xg = (min(x_gridding):bostick_grid(1):max(x_gridding));               
% yg = (min(y_gridding):bostick_grid(2):max(y_gridding));
xg = (maplim(1):bostick_grid(1):maplim(2));               
yg = (maplim(3):bostick_grid(2):maplim(4));
[X,Y] = meshgrid(xg,yg);  

%==========================================================================
for i=1:size(bostick_pseudo_av,1)
    bostick_rho_gridding = log10(bostick_pseudo_av(i,:));
    x_gridding=x_sort;     y_gridding=y_sort;
    % delete the data corresponding to NAN in bostick_pseudo
    x_gridding(isnan(bostick_rho_gridding))=[];
    y_gridding(isnan(bostick_rho_gridding))=[];
    bostick_rho_gridding(isnan(bostick_rho_gridding))=[];
    if length(x_gridding)<=5
        bostick_grid_av(:,:,i) = nan(size(X));
    else
%         bostick_grid_av(:,:,i) = griddata(x_gridding,y_gridding,bostick_rho_gridding,X,Y,'cubic');
        F = scatteredInterpolant(x_gridding',y_gridding',bostick_rho_gridding','natural','linear');
        % You can use natural or nearest as the interpolation method, and
        % nearest is closer to the real data but rougher than the natural
        % method.
        bostick_grid_av(:,:,i) = F (X,Y);
    end
end

%==========================================================================
% Plot the Bostick results at depths defined in bostick_depth

for i=1:length(bostick_depth)
    dtemp=abs(d_interp-log10(1000*bostick_depth(i)));
    nearest=min(dtemp);
    index=find(dtemp==nearest);
    
    figure(fig_number);
    [h]=m_pcolor(xg,yg,bostick_grid_av(:,:,index));hold on;
    m_grid('tickdir','out','yaxislocation','left','xlabeldir','end','ticklen',.02);
    m_plot(x,y,'k.');    
%   m_plot(xg,yg,'r.')
    plot_geo_bound_file_m_map(geo_bound_file);
    
    caxis([log10(rlim(1)),log10(rlim(2))]);    
    colormap(flipud(jet));    colorbar vertical;
    shading flat;    hold off;
    % xlabel('longitude'); ylabel('latitude');  
    % calculate the real depth closest to the given depth
    depth=(10^d_interp(index))/1000;
    title(['Average Bostick resistivity (ohm.m); z = ',num2str(depth,'%5.1f'),'  km;', ...
           ' dx =',num2str(bostick_grid(1),'%3.0f'),' km; dy =',num2str(bostick_grid(2),'%3.0f'),' km'])

    % output the figure
    eval(['print -djpeg100 bostick\bostick_gridding1_layer',...
        num2str(i,'%03.0f'),'_depth_',num2str(round(depth),'%5.1f\n'),'_km.jpg']);
    
 
    % output Bostick resistivity to a text file (use in Surfer or ArcGIS)
    str=['bostick\bostick_av_layer_',num2str(i,'%03.0f'),'_depth_',num2str(round(depth),'%5.1f\n'),'_km.dat'];
    fid2=fopen(str,'w+');
    for ix=1:length(xg)
      for iy=1:length(yg)
        fprintf(fid2,'%9.3f %9.3f %9.3f %9.3f \n',[xg(ix),yg(iy),depth,bostick_grid_av(iy,ix,index)]);
      end
    end
    fclose(fid2);
    
    pause (1)
    close(fig_number)

end 
end