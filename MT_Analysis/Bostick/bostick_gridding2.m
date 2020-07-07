function [p1] = bostick_gridding2(fig_number,Z,Zvar,T,azimuth,rlim,bostick_grid,bostick_int_z,bostick_depth,stn_name,lat,long,elev,maplim,ns,profile_dir,geo_bound_file)

% 
% Plots a map of the Bostick resistivity at a set of depths for a grid
% of MT stations.  This is a advanced version, calculating the Bostick
% resistivity from the Berdichevsky aveage and direction of rho_max

% bostick_gridding1.m is simpler (faster) version

rad = 180./pi;      p1 = 1;   prof_azim = mean(azimuth);

mkdir bostick

% Process co-ords
x = long; y=lat;
% x_mean = mean(x);  y_mean = mean(y); % km
% x = cos(y_mean/rad)*111*(x-x_mean);
% y = 111*(y-y_mean); % km


c = cosd(prof_azim);      s = sind(prof_azim);
R = [ c, -s ; s, c];

for is=1:ns
    loc = [x(is),y(is)];
    loc = R*loc';
    x_rot(is) = loc(1);y_rot(is)=loc(2);
end
  
% Add loop to put stations in order on monotonically increasing x_rot
[x_sort,index] = sort(x_rot,'ascend');
  


     
%==========================================================================
% calculate the bostick from the average
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
logdepth=log10(d_av);
d_interp=min(min(logdepth)):bostick_int_z:max(max(logdepth));% Define depths onto which we interpolate

% Set up the grid  (in degree) for lateral interpolation 
% x_gridding=x_sort;     y_gridding=y_sort;
% xg = (min(x_gridding):bostick_grid(1):max(x_gridding));               
% yg = (min(y_gridding):bostick_grid(2):max(y_gridding));
xg = (maplim(1):bostick_grid(1):maplim(2));               
yg = (maplim(3):bostick_grid(2):maplim(4));
[X,Y] = meshgrid(xg,yg);  

%==========================================================================
% Calculate the maximum resistivity corresponding to each layer in d_interp
for is=1:ns
    i_sort = index(is);  
    Z1 = Z(:,:,:,i_sort);    Z1var = Zvar(:,:,:,i_sort);  
    [~,~,~,~,~,bostick_pseudo_max(:,is),~,~,bostick_max_dir(:,is),~] = bostick2(Z1,Z1var,T,d_interp);
end

%==========================================================================
for i=1:size(bostick_pseudo_max,1)
    bostick_rho_gridding = log10(bostick_pseudo_max(i,:));
    x_gridding=x_sort;     y_gridding=y_sort;
    % delete the data corresponding to NAN in bostick_pseudo
    x_gridding(isnan(bostick_rho_gridding))=[];
    y_gridding(isnan(bostick_rho_gridding))=[];
    bostick_rho_gridding(isnan(bostick_rho_gridding))=[];
    if length(x_gridding)<=5
        bostick_grid_max(:,:,i) = nan(size(X));
    else
%         bostick_grid_max(:,:,i) = griddata(x_gridding,y_gridding,bostick_rho_gridding,X,Y,'cubic');
%         F = TriScatteredInterp(x_gridding',y_gridding',bostick_rho_gridding','natural');
        F = scatteredInterpolant(x_gridding',y_gridding',bostick_rho_gridding','natural','linear');
        % You can use natural or nearest as the interpolation method, and
        % nearest is closer to the real data but rougher than the natural
        % method.
        bostick_grid_max(:,:,i) = F (X,Y);
    end
end

%==========================================================================
% Plot the Bostick results at depths defined in bostick_depth

for i=1:length(bostick_depth)
    dtemp=abs(d_interp-log10(1000*bostick_depth(i)));
    nearest=min(dtemp);
    index=find(dtemp==nearest);
    
    figure(fig_number);
    m_grid('tickdir','out','yaxislocation','left','xlabeldir','end','ticklen',.02);
    hold on;  
    m_pcolor(xg,yg,bostick_grid_max(:,:,index)); 
    m_plot(x,y,'k.'); 
    plot_geo_bound_file_m_map(geo_bound_file);
    
    caxis([log10(rlim(1)),log10(rlim(2))]);
    colormap(flipud(jet));    colorbar vertical;
    shading flat;   
%     axis([min(xg),max(xg),min(yg),max(yg)]);  
    % calculate the real depth closest to the given depth
    depth=(10^d_interp(index))/1000;
    title(['Max Bostick resistivity (ohm.m); z = ',num2str(depth,'%4.0f'),'  km;', ...
           ' dx =',num2str(bostick_grid(1),'%3.0f'),' km; dy =',num2str(bostick_grid(2),'%3.0f'),' km'])
%     xlabel('longitude'); ylabel('latitude'); 
    
    % plot the direction of the rho_max
    max_dir_gridding = (bostick_max_dir(index,:));
    x_gridding=x_sort;     y_gridding=y_sort;
    % delete the data corresponding to NAN in bostick_pseudo
    x_gridding(isnan(max_dir_gridding))=[];
    y_gridding(isnan(max_dir_gridding))=[];
    max_dir_gridding(isnan(max_dir_gridding))=[];
    
    e_comp=cosd(90-max_dir_gridding);  % east componant of the arrows
    n_comp=sind(90-max_dir_gridding);  % north componant of the arrows    
   
    m_quiver(x_gridding,y_gridding,e_comp,n_comp,0.3,'k');
    m_quiver(x_gridding,y_gridding,-1*e_comp,-1*n_comp,0.3,'k');
    hold off;

    % output the figure
    eval(['print -djpeg100 bostick\bostick_gridding2_depth_',...
        num2str(round(depth),'%3.3i\n'),'_km.jpg']);
    
    pause(2)
    close(fig_number)

end 
end