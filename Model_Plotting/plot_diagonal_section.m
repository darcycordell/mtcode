function [x,y,z,rho] = plot_diagonal_section(m,d)
% Function to plot diagonal cross-sections through a 3D resistivity model
% User chooses the section to plot on the map slice or by inputing lat-long
% points
%
% Usage: [x,y,z,rho] = plot_diagonal_section(m,d)
%   OR   [x,y,z,rho] = plot_diagonal_section(m)
%
% Inputs:
% "m" is the model structure
%
% "d" is OPTIONAL to data structure to plot site locations on sections
%
% Outputs:
% x,y,z are the NS, EW and depth locations of the section
% rho is the matrix of resistivity values associated with the section
%
% This function uses the "line" method where it calculates line segment
% distances through each cell along the diagonal trace and plots the
% resistivity values for each of those line segments. This results in a
% pcolor plot with different sized cells because the diagonal trace may
% only just clip the edges of some cells. See
% plot_diagonal_section_smooth.m for alternative method.
%


u = user_defaults;
% [L] = load_geoboundary_file_list; % plot_geoboundaries_diagonal_section
% calls this
close all

%If d structure does not exist then set the variable to all NaN
if ~exist('d','var')
    [d] = make_nan_data;
end

        
%------------------Plot a Slice to Choose Diagonal Section------------------
id = round(m.nz/2); %index to plot a slice (arbitrary)
figure(1);
plot_slice(m,id,d); hold on;
title(['Depth = ',num2str(m.cz(id)/1000),' km b.s.l.']);

%Option to pick points by clicking or manually entering lat-long points
%Clicking is more common, but for repeatability, if you want to specify two
%points which you had previous clicked, then you can enter them manually as
%well.
main_menu=menu('Define diagonal slice by:','Clicking','Enter X/Y Coordinates Manually','Load Text File With X/Y Endpoints','Enter Lat/Long Coordinates Manually','Load Text File With Lat/Lon Endpoints');

if main_menu == 1 %Option to Click Points

    disp('Choosing:  Diagonal slice')
    disp('Choose points to define a line or fence section, right click to exit')
    i = 1; xp = [0,0]; yp = [0,0];
    while i<=2
        [yp(i), xp(i)]=ginput(1);  %graphically choose the slice to be plotted
        plot(yp(i),xp(i),'r*');
        i = i+1;
    end
    
    if isstruct(d) %If data exists, then plot station locations
        [lon,lat]=utm2geo(yp*1000+500000, xp*1000,d.origin(2), d.origin(1));
    else
        lon = yp;
        lat = xp;
    end
    
    plot(yp,xp,'-r'); %Plot the slice in map view that will be plotted

    %Code displays the lpoints pairs which were clicked for repeatability.
    %The next time, if you want to have the identical slice, just manually
    %enter the coordinates rather than clicking.
    disp('You chose the points:')

    for i = 1:length(lon)
   
        disp(['(',num2str(lat(i)),', ',num2str(lon(i)),')'])
    end
    
    %These are the chosen points in meters
    xp = xp*1000;
    yp = yp*1000;
    
elseif main_menu == 2 % enter X/Y end points manually
    
    prompt={sprintf('Enter in km\n\nNS 1'),'EW 1','NS 2','EW 2'}; 
    dlg_title='Kilometers';
    def={'110','-405','-10.789','63.8'}; %default is office to Tim's
    num_lines=1;
    dinp = inputdlg(prompt,dlg_title,num_lines,def);

    xp=[str2double(dinp{1});str2double(dinp{3})]; 
    yp=[str2double(dinp{2});str2double(dinp{4})];

    %These are the chosen points
    xp = xp*1000;
    yp = yp*1000;
    
    plot(yp/1000,xp/1000,'-*r'); %Plot the diagonal trace in map view that will be plotted
    
elseif main_menu == 3 % load X/Y endpoints from a text file
    
    curdir = pwd;
    [profile_points_file, filepath]=uigetfile({'*.txt'},'Choose text file which contains profile endpoints'); 
    if profile_points_file == 0 
        return; 
    end
    cd(filepath)

    profile_loaded=load(profile_points_file);
    fclose all;
    cd(curdir)
    yp = profile_loaded(:,1); % transpose just to keep dimensions consistent with clicking on map
    xp = profile_loaded(:,2); 
    
    xp = xp*1000;
    yp = yp*1000;
    
    plot(yp/1000,xp/1000,'-*r'); %Plot the diagonal trace in map view that will be plotted
    
elseif main_menu == 4 %Option to enter lat/long pair manually
    
    if unique(d.zrot) ~= 0 % data are rotated, assume model correpsonds to the rotated data set
        error('Warning: profile only working in model coordinates for rotated models. Use options 1, 2, or 3.')
    end
    
    if ~strcmp(d.site,'None') %If a data structure was supplied, then get the points in latitude and longitude
        prompt={'Latitude 1','Longitude 1','Latitude 2','Longitude 2'}; 
        dlg_title='Decimal Degrees';
        def={num2str(d.origin(1)),num2str(d.origin(2)),num2str(d.origin(1)),num2str(d.origin(2))}; %default is office to Tim's
        num_lines=1;
        dinp = inputdlg(prompt,dlg_title,num_lines,def);

        lat=[str2double(dinp{1});str2double(dinp{3})]; 
        lon=[str2double(dinp{2});str2double(dinp{4})];

        %Convert latitude and longitudes to x-y model coordinates
        [yp,xp] = geo2utm(lon,lat,d.origin(2), d.origin(1));
        
        yp = yp -500000;
        
    else %If no data structure was supplied, then get the points in x-y model coordinates
        
        prompt={'NS 1','EW 1','NS 2','EW 2'}; 
        dlg_title='Kilometers';
        def={'110','-405','-10.789','63.8'}; %default is office to Tim's
        num_lines=1;
        dinp = inputdlg(prompt,dlg_title,num_lines,def);

        xp=[str2double(dinp{1});str2double(dinp{3})]; 
        yp=[str2double(dinp{2});str2double(dinp{4})];
        
        %These are the chosen points
        xp = xp*1000;
        yp = yp*1000;
    
    end

    plot(yp/1000,xp/1000,'-*r'); %Plot the diagonal trace in map view that will be plotted
    
elseif main_menu == 5 %Option to load text file with profile endpoints
        
    if unique(d.zrot) ~= 0 % data are rotated, assume model correpsonds to the rotated data set
        error('Warning: profile only working in model coordinates for rotated models. Use options 1, 2, or 3.')
    end
    
    curdir = pwd;
    [profile_points_file, filepath]=uigetfile({'*.txt'},'Choose text file which contains profile endpoints'); 
    if profile_points_file == 0 
        return; 
    end
    cd(filepath)

    profile_loaded=load(profile_points_file);
    fclose all;
    cd(curdir)
    lon = profile_loaded(:,1); % transpose just to keep dimensions consistent with clicking on map
    lat = profile_loaded(:,2); 

    [yp,xp] = geo2utm(lon,lat,d.origin(2), d.origin(1));
    yp = yp -500000;    
    
    plot(yp/1000,xp/1000,'-*r'); %Plot the slice in map view that will be plotted   
    
else
    
    return
    
end

%Sort the points so that profile is always plotted the same way regardless
%of the order you pick the points.
[xp, id] = sort(xp,'descend');
yp = yp(id);
%%
if max(xp)>max(m.cx) || min(xp)<min(m.cx) || max(yp)>max(m.cy) || min(yp)<min(m.cy)
    disp('Error: Your lat/long ranges are outside the model space')
else

    %Step 1: Find the slope of the line and intercept of the line between the
    %two points

    s=(xp(2)-xp(1))/(yp(2)-yp(1)); %slope
    b=xp(1)-s*yp(1); %intercept

    if isinf(s)
        disp('The diagonal selection does not work for NS-profiles. Choose to plot NS section in S3 instead');
        return;
    end


    %Step 2: Find the intersections between the line and the grid edges in NS
    %and EW directions

    %EW first (xint: dealing with known east-west points)
    int=[];
    for k=1:m.ny
        xint=s*m.cy(k)+b; %slope formula to find intersections with known x
        if m.cy(k)>=min(yp) && m.cy(k)<=max(yp)
            int=[int; [m.cy(k) xint]];
        end
    end

    %NS second (xint: dealing with known north-south points)
    for k=1:m.nx
        yint=(m.cx(k)-b)/s; %slope formula to find intersections with known y
        if m.cx(k)>=min(xp) && m.cx(k)<=max(xp)
            int=[int; [yint m.cx(k)]];
        end
    end

    %add beginning and end points to intercept matrix and sort
    int=[int;[yp(1) xp(1)];[yp(2) xp(2)]];
    
    %Sort the intercepts so that profile is always plotted the same way regardless
    %of the slope of the line. It is plotted such that the WEST side of
    %the profie is on the LEFT and the EAST side of the profile is on the
    %RIGHT. Similarly, the NORTH side of the profile is on the LEFT and the
    %SOUTH side of the profile is on the RIGHT.
    [int, idx] =sortrows(unique(int,'rows'),2); %remove any duplicates
    if s<0 %For sections with negative slopes (i.e. NW-SE)
        int = int(idx,:);
    else
        if s>=1 %For sections with positive slopes >45 degrees (i.e. steep NE-SW)
            int = int(flipud(idx),:); 
        else %For sections with positive slopes <45 degrees (i.e. shallow NE-SW)
            int = int(idx,:);
        end
    end
    
    %Option to plot the intercepts on the map (for debugging purposes)
%     plot(int(:,1)/1000,int(:,2)/1000,'*w')
%     plot(int(1,1)/1000,int(1,2)/1000,'*k','MarkerSize',12)

    %Step 3: Find the distance and midpoint of each line segment through each grid cell

    %Get lengths and midpoints of line segments
    L=zeros(length(int(:,1))-1,1);L_mid=zeros(length(L),2); L_d = zeros(size(L));
    for k=1:length(int)-1 %loops over all grid intersections
        L_set=sqrt((int(k+1,1)-int(k,1))^2+(int(k+1,2)-int(k,2))^2);
        L(k)=L_set;
        L_mid_set=(int(k+1,:)+int(k,:))./2;
        L_mid(k,:)=L_mid_set;

        L_d(k) = L_set;

    end

    dist=cumsum(L_d);

    %Step 4: Find which cell the midpoint belongs in
    x_grid = zeros(length(L),1); y_grid = zeros(length(L),1);
    for j=1:length(L)
        y_grid(j) = nearestpoint(L_mid(j,1),m.cy);
        x_grid(j) = nearestpoint(L_mid(j,2),m.cx);
    end
    
    %Option to plot the grid cells on the map (for debugging purposes)
%     plot(m.cy(y_grid)/1000,m.cx(x_grid)/1000,'xk')
%     plot(m.cy(y_grid(1))/1000,m.cx(x_grid(1))/1000,'xr','MarkerSize',12)
%%
    res_plot=zeros(m.nz,length(L));
    for j=1:length(L)
        res_plot(:,j) = m.A(x_grid(j),y_grid(j),:); %pulls slice model cells from A matrix
    end
        
    if strcmp(u.diagonal_section_mode,'interp')
        ind = find(isnan(res_plot));
        res_plot = inpaint_nans(res_plot);
        px = cumsum(ones(length(L),1)*sum(L)/(length(L)));
        [D,Zorig] = meshgrid(dist,m.cz);
        [P,Z] = meshgrid(px,m.cz);
        res_plot =interp2(D,Zorig,res_plot,P,Z);
        res_plot(ind) = NaN;
        dist = px;
    end
    %%
    if ~strcmp(d.site,'None')
        %Project sites within the tolerance distance from the profile onto the
        %profile
        k=0; id = []; d_elev = [];
        for i = 1:length(d.idx)
            

                %Find the minimum distance between the site and the profile cell
                %centers. This is the cell center which the site will be
                %projected to
                [d_site, ind] = min(sqrt((m.cx(d.idx(i))-L_mid(:,2)).^2+(m.cy(d.idy(i))-L_mid(:,1)).^2));

                if d_site<u.tol*1000 %if station is less than tolerance (km) from profile then include it on the profile plot
                    k=k+1;
                    id(k) = i; %Projected station index
                    
                    %d_elev includes the model indices of the projected
                    %site and the distance the site is along the profile
                    d_elev(k,:) = [x_grid(ind) y_grid(ind) dist(ind)];
                end
            
        end
        
        %Plot which sites will be projected to the profile
        figure(1)
        plot(d.y(id)/1000,d.x(id)/1000,'r.','markersize',14);   
        %Option to plot the cell centers which the stations are projected to
        %plot(m.cy(d_elev(:,2))/1000,m.cx(d_elev(:,1))/1000,'xr')
  
    end
    
    zind = 1:nearestpoint(u.zmax*1000,m.cz);

    %Plot section
    figure(2);
   
    plot_cross_section(dist/1000,m.cz(zind)/1000,res_plot(zind,:)'); % transpose to negate transpose in this function...
    
    xlabel(['Kilometers from (',num2str(xp(1)/1000),' NS, ',num2str(yp(1)/1000),' EW) to (',num2str(xp(2)/1000),' NS, ',num2str(yp(2)/1000),' EW)']);
    ylabel('Depth Below Sealevel (km)')
    title('Diagonal Slice Through Model');

    %%
%     seisy = dist/1000;
%     seisz = m.cz(zind)/1000;
%     seisA = res_plot(zind,:);
%     save seis_2d_for_interpolation_v2 seisy seisz seisA
    
    %Plot geoboundary file on diagonal section (not yet implemented)
    plot_geoboundaries_diagonal_section(d,L_mid(:,1),L_mid(:,2))    
    
    %Plot projected station locations (projected onto the profile and
    %projected onto the elevation surface of the profile)
    if ~strcmp(d.site,'None')
        if ~isempty(d_elev)
            figure(2);
            elev = m.Z(sub2ind(size(m.Z),d_elev(:,1),d_elev(:,2)));
            plot(d_elev(:,3)/1000,elev/1000 -u.zoff,'vk','MarkerFaceColor','k'); hold on
        end
    end
    

    print_figure(['diagonal_slices_',m.niter],['slice_from_',num2str(xp(1)/1000),'_',num2str(yp(1)/1000),'_to_',num2str(xp(2)/1000),'_',num2str(yp(2)/1000)])
    
    figure(1)
    print_figure(['diagonal_slices_',m.niter],['slice_from_',num2str(xp(1)/1000),'_',num2str(yp(1)/1000),'_to_',num2str(xp(2)/1000),'_',num2str(yp(2)/1000),'_MAP']) 


    %Output x,y,z,rho which is the x,y,z coordinates and corresponding rho
    %values of the plotted section
    x = L_mid(:,2)/1000;
    y = L_mid(:,1)/1000;
    z = m.cz(zind)/1000;
    rho = res_plot(zind,:)';

end
end %End function