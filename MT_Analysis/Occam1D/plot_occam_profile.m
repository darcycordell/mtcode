clear all; close all

%-------------------------------------------------------------------------
%INPUTS:
% After running OCCAM1D for each EDI file, you have the option to save the
% inversion output as a *.mat file.
%
% List of *.mat files with OCCAM inversion results:
files = {'ARG001.mat','ARG002.mat','ARG003.mat','ARG004.mat'};

reslims = [1 1000]; %Resistivity limits for plotting
depthlims = [0 5000]; %Depth limits for plotting (in meters)

%-------------------------------------------------------------------------


%Load OCCAM1D results
for i = 1:length(files)
    inv(i) = load(files{i});
end

%Plot map of site locations:
lon = [inv(:).lon]; lat = [inv(:).lat];
plot(lon,lat,'*k'); hold on
xlabel('Longitude'); ylabel('Latitude')

%Get "profile" through the stations
L = polyfit(lon,lat,1);
x = min(lon):0.001:max(lon);
y = polyval(L,x);
plot(x,y,'--k')

id = dsearchn([x;y]',[lon; lat]'); %Profile indices where stations are
%plot(x(id),y(id),'*r') %For debugging to check

%Get distance each site is along the profile
dist = [];
for i = 1:length(inv)
    dist(i) = haversine(y(1),x(1),y(id(i)),x(id(i)));
end

%Plot figure with the profile distance as the x axis
set_figure_size(2); h = plot([dist(1)-max(diff(dist))*0.1,dist(end)+max(diff(dist))*0.1],[0 0],'-w'); hold on
pos = get(gca, 'Position');
ax1 = gca;
set(ax1,'YTick',[])
xlabel('Distance Along Profile (km)')
fig1ylim = ylim;
fig1xlim = xlim;

%Plot each 1-D model as an inset figure at the correct x position.
for i = 1:length(inv)
    plot(ax1,dist(i),0.8,'vk','MarkerFaceColor','k','MarkerSize',12)
    ax(i) = axes('Position',[(dist(i))/diff(fig1xlim)*pos(3)+pos(1)+(max(diff(dist))*0.1/diff(fig1xlim)), ...
        pos(2)+0.1*pos(4), ...
        1/diff(fig1ylim) * pos(3)/length(inv), ...
        pos(4)-2*0.08*pos(4)]);
    
    %ax3 = axes('Position',[pos(1) pos(2) pos(3) pos(4)])
    stairs(inv(i).model,inv(i).depth/1000,'-k','LineWidth',1);
    set(ax(i),'XScale','log')
    xlabel('Resistivity (\Omega m)')
    ylabel('Depth (km)')
    axis([reslims depthlims/1000])
    axis ij
    grid on
    
end

