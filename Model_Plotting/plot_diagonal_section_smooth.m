function [x,y,z,rho] = plot_diagonal_section_smooth(m,d)
% Function to plot diagonal cross-sections through a 3D resistivity model
% User chooses the section to plot on the map slice or by inputing lat-long
% points
%
% Usage: [x,y,z,rho] = plot_diagonal_section_smooth(m,d)
%   OR   [x,y,z,rho] = plot_diagonal_section_smooth(m)
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
% This function uses the "smooth" method where it calculates the diagonal
% section by smoothing over shared cells so that the plotted section has 
% uniform cell sizes. For sections which are near 45 degree slope, this 
% method works very well. But for sections near 0 or 90 degrees, it will 
% not look good. See plot_diagonal_section.m for alternative method.

u = user_defaults;
[L] = load_geoboundary_file_list;
close all

%If d structure does not exist then set the variable to all NaN
if ~exist('d','var')
    [d] = make_nan_data;
end
        
%------------------Plot a Slice to Choose Diagonal Slice------------------
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
    
    xp = xp*1000;
    yp = yp*1000;
    
elseif main_menu == 2 % enter X/Y end points manually
    
    prompt={sprintf('Enter in km\n\nNS 1'),'EW 1','NS 2','EW 2'}; 
    dlg_title='Kilometers';
    def={'0','-5','5','10'};
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

        [yp,xp] = geo2utm(lon,lat,d.origin(2), d.origin(1));
        
        yp = yp -500000;
        
    else
        
        prompt={'NS 1','EW 1','NS 2','EW 2'}; 
        dlg_title='Kilometers';
        def={'110','0','0','0'}; %default is office to Tim's
        num_lines=1;
        dinp = inputdlg(prompt,dlg_title,num_lines,def);

        xp=[str2double(dinp{1});str2double(dinp{3})]; 
        yp=[str2double(dinp{2});str2double(dinp{4})];
        
        %These are the chosen points in meters
        xp = xp*1000;
        yp = yp*1000;

    end

    plot(yp/1000,xp/1000,'-*r'); %Plot the slice in map view that will be plotted
 
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

if max(xp)>max(m.cx) || min(xp)<min(m.cx) || max(yp)>max(m.cy) || min(yp)<min(m.cy)
    disp('Note: Your lat/long ranges are outside the model space')
else

    %Step 1: Find the slope of the line and intercept of the line between the
    %two points

    s=(xp(2)-xp(1))/(yp(2)-yp(1)); %slope
    b=xp(1)-s*yp(1); %intercept

    if isinf(s)
        disp('The smooth diagonal selection does not work for NS-profiles. In S3, choose "NS slice"');
        return;
    end
    
    if s==0
        disp('The smooth diagonal selection does not work for EW profiles. In S3, choose "EW Slice"')
        return
    end
    
    


    %Step 2: Find the intersections between the line and the grid edges in NS
    %and EW directions

    %EW first (xint: dealing with known east-west points)
    int=[];
    for k=1:m.ny+1
        xint=s*m.y(k)+b; %slope formula to find intersections with known x
        if m.y(k)>=min(yp) && m.y(k)<=max(yp) && ~isnan(xint)
            int=[int; [m.y(k) xint]];
        end
    end

    %NS second (xint: dealing with known north-south points)
    for k=1:m.nx+1
        yint=(m.x(k)-b)/s; %slope formula to find intersections with known y
        if m.x(k)>=min(xp) && m.x(k)<=max(xp) && ~isnan(yint)
            int=[int; [yint m.x(k)]];
        end
    end

    %add beginning and end points to intercept matrix and sort
    int=[int;[yp(1) xp(1)];[yp(2) xp(2)]];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Find/average the resistivity blocks along the profile 
ns_ind = nearestpoint(int(:,2),m.cx);
ew_ind = nearestpoint(int(:,1),m.cy);

if length(unique(ns_ind))>length(unique(ew_ind))  %If the section is >45 degrees
    
    [ns_ind, indn] = unique(ns_ind);
    ew_ind = ew_ind(indn);
    
    %Plot the indices to be averaged (for debugging purposes)
    %plot(m.cy(ew_ind)/1000,m.cx(ns_ind)/1000,'xr')
    
    [ind_lean, ia] =unique(ew_ind);
    ew_ind_A = ew_ind(ia);
    ns_ind_A = ns_ind(ia);
    
    res_plot=nan(length(m.cz),length(ind_lean));
    ew_locations = zeros(length(ind_lean),1);
    ns_locations = zeros(length(ind_lean),1);
    for i=1:length(ind_lean)

        ind_temp=find(ew_ind==ind_lean(i));

        if length(ind_temp)>1

            res_sum=zeros(1,1,m.nz);
            for j=1:length(ind_temp)
                res_sum=res_sum+m.A(ns_ind(ind_temp(j)),ew_ind(ind_temp(j)),:);
            end

            res_plot(:,i)=squeeze(res_sum)./(length(ind_temp));
            ew_locations(i)=m.y(ew_ind(ind_temp(1)));
            ns_locations(i)=mean(m.x(ns_ind(ind_temp)));
            

        else

            res_plot(:,i)=squeeze(m.A(ns_ind(ind_temp),ew_ind(ind_temp),:));
            ew_locations(i)=m.y(ew_ind(ind_temp));
            ns_locations(i)=m.x(ns_ind(ind_temp));
           
        end       
    end   
    
else %If the section is less than 45 degrees slope
 
    [ew_ind, inde] = unique(ew_ind);
    ns_ind = ns_ind(inde);
    
    %Plot the locations of the cells to be used (for debugging purposes)
    %plot(m.cy(ew_ind)/1000,m.cx(ns_ind)/1000,'xw')
    
    [ind_lean, ia] =unique(ns_ind);
    ns_ind_A = ns_ind(ia);
    ew_ind_A = ew_ind(ia);
    
    res_plot=nan(length(m.cz),length(ind_lean));
    
    ew_locations = zeros(length(ind_lean),1);
    ns_locations = zeros(length(ind_lean),1);
    for i=1:length(ind_lean)

        ind_temp=find(ns_ind==ind_lean(i));
        
        if length(ind_temp)>1

            res_sum=zeros(1,1,m.nz);

            for j=1:length (ind_temp)
                res_sum=res_sum+m.A(ns_ind(ind_temp(j)),ew_ind(ind_temp(j)),:);
            end

            res_plot(:,i)=squeeze(res_sum)./(length(ind_temp));
            ns_locations(i)=m.cx(ns_ind(ind_temp(1)));
            ew_locations(i)=mean(m.cy(ew_ind(ind_temp)));
        else
            res_plot(:,i)=squeeze(m.A(ns_ind(ind_temp),ew_ind(ind_temp),:));
            ns_locations(i)=m.cx(ns_ind(ind_temp));
            ew_locations(i)=m.cy(ew_ind(ind_temp));
        end        
    end
end

%Plot the cell locations (for debugging purposes)
%plot(ew_locations/1000,ns_locations/1000,'*k')


%find the distance of cells in the profile from one end 
length_segments = zeros(length(ns_locations)-1,1);
for i= 1: (length(ns_locations)-1 )
    length_segments(i)=sqrt((ns_locations(i+1)-ns_locations(i))^2+...
        (ew_locations(i+1)-ew_locations(i))^2);
end

dist = [0; cumsum(length_segments)];

%Sort the distances so that profile is always plotted the same way regardless
%of the slope of the line. It is plotted such that the WEST side of
%the profie is on the LEFT and the EAST side of the profile is on the
%RIGHT. Similarly, the NORTH side of the profile is on the LEFT and the
%SOUTH side of the profile is on the RIGHT.
if s>1
    dist = flipud(dist);
elseif s<0 && s>-1
    dist = flipud(dist);
end


if ~strcmp(d.site,'None') %If data exists then determine which sites should be put on the profile and where
        
    %Project sites within the tolerance distance from the profile onto the
    %profile

    k=0; id = []; d_elev = [];
    for i = 1:length(d.idx)
        for j = 1:length(ind_lean)

            %Find the distance between the sites and the profile
            d_site = sqrt((m.cx(d.idx(i))-ns_locations(j))^2+(m.cy(d.idy(i))-ew_locations(j))^2);

            if d_site<u.tol*1000 %if station is less than tolerance from profile then include it on the profile plot
                k=k+1;
                id(k) = i;
                d_elev(k,:) = [ns_ind_A(j) ew_ind_A(j)];
            end
        end
    end

    [id,ia] = unique(id); %find all unique sites
    d_elev = d_elev(ia,:);

    figure(1)
    plot(d.y(id)/1000,d.x(id)/1000,'r.','markersize',14);

    %Organize the site distances so that they are plotted correctly
    %from west to east (or north to south)
    if s<0 %For NW-SE profiles
        d_site = sqrt((xp(1)-d.x(id)).^2+(yp(1)-d.y(id)).^2);
    else
        if s<1 %For steep NE-SW profiles >45 degree slope
            d_site = sqrt((xp(2)-d.x(id)).^2+(yp(2)-d.y(id)).^2);
        else %For shallow NE-SW profile <45 degree slope
            d_site = sqrt((xp(1)-d.x(id)).^2+(yp(1)-d.y(id)).^2);
        end

    end

end

zind = 1:nearestpoint(u.zmax*1000,m.cz);

%Plot section
figure(2);
plot_cross_section(dist/1000,m.cz(zind)/1000,log10(res_plot(zind,:))); % transpose to negate transpose in this function...

xlabel(['Kilometers from (',num2str(xp(1)/1000),' NS, ',num2str(yp(1)/1000),' EW) to (',num2str(xp(2)/1000),' NS, ',num2str(yp(2)/1000),' EW)']);
ylabel('Depth Below Sealevel (km)')
title('Diagonal Slice Through Model');

if strcmp(u.gridlines,'off')
    shading flat
end

%Plot geoboundary file on diagonal section (not yet implemented properly)
plot_geoboundaries_diagonal_section(d,ew_locations,ns_locations)


if ~strcmp(d.site,'None') %If data exists plot them
    if ~isempty(d_elev)
        figure(2);
        elev = m.Z(sub2ind(size(m.Z),d_elev(:,1),d_elev(:,2)));
        plot(d_site/1000,elev/1000 -u.zoff,'vk','MarkerFaceColor','k'); hold on
    end
end


print_figure(['diagonal_slices_',m.niter],['smooth_slice_from_',num2str(xp(1)/1000),'_',num2str(yp(1)/1000),'_to_',num2str(xp(2)/1000),'_',num2str(yp(2)/1000)]) 


figure(1)
print_figure(['diagonal_slices_',m.niter],['smooth_slice_from_',num2str(xp(1)/1000),'_',num2str(yp(1)/1000),'_to_',num2str(xp(2)/1000),'_',num2str(yp(2)/1000),'_MAP'])

%Output x,y,z,rho which is the x,y,z coordinates and corresponding rho
%values of the plotted section
x = ns_locations/1000;
y = ew_locations/1000;
z = m.cz(zind)/1000;
rho = res_plot(zind,:)';


  
end

end