function plot_skew_pseudo(d)
%
% Function which plots pseudo section of Swift skew and Bahr skew
%
% Usage: plot_skew_pseudo(d)
%
% "d" is an MT data structure
%
u = user_defaults;

%Pull a profile to plot the pseudo-section. "sidx" is stations indices
%included on profile and "rot_ang" is the direction of the normal to the
%profile direction
[sidx,midpoint,azimuth] = get_pseudo_section_indices(d);

rad = 180./pi;      

% Process cooordinates
lon_mean = midpoint(1);  lat_mean = midpoint(2); % km
x = cos(lat_mean/rad)*111*(d.loc(sidx,2)-lon_mean);y = 111*(d.loc(sidx,1)-lat_mean); % km
N = length(sidx);


c = cosd(azimuth);      s = sind(azimuth);     R = [ c, -s ; s, c];

%Rotate station coordinates to get distance along profile
x_rot = zeros(N,1); y_rot = zeros(size(x_rot));
for is=1:N
    loc = [x(is),y(is)];    loc = R*loc';
    x_rot(is) = loc(1);     y_rot(is)=loc(2);
end

% Add loop to put stations in order on monotonically increasing x_rot
[x_sort,index] = sort(x_rot,'ascend');
index = sidx(index);

if ~u.rotate_data_to_azimuth
    dr=d;
else
    [dr]=rotate_d(d,azimuth); %Rotate data   
end

for is = 1:N
    [dim] = calc_dim_parameters(dr.Z(:,:,index(is)),azimuth);
    
    swift_skew(:,is) = dim.swift_skew;
    bahr_skew(:,is) = dim.bahr_skew;
    
end
%%
%========================================================================
%Plot swift skew
set_figure_size(1);
colormap(flipud(u.cmap));
subplot(2,1,1)
[hp,vp,C] = fix_pcolor(x_sort,log10(d.T),swift_skew);
pcolor(hp,vp,C)
hold on; axis ij
shading 'flat'
plot(x_sort,log10(d.T(1))-0.25,'kv','MarkerFaceColor','k')
if u.station_names
    text(x_sort,(x_sort.*0 )+ min(vp)-0.5,d.site(index),'rotation',u.station_names_angle,'interpreter','none','fontsize',8);
else % station names and title overlap. if station name shown, quanitity still shown in colorbar label
    title('Swift Skew')
end

ylabel ('Log10 Period (s)')
axis([hp(1),hp(end), min(vp)-0.5 ,max(vp)])
caxis([0 0.4])
hcb = colorbar;
hcb.Label.String = 'Swift Skew';
set(gca,'Layer','top')

%Plot Bahr skew
subplot(2,1,2)
[hp,vp,C] = fix_pcolor(x_sort,log10(d.T),bahr_skew);
pcolor(hp,vp,C)
hold on; axis ij
shading 'flat'
plot(x_sort,log10(d.T(1))-0.25,'kv','MarkerFaceColor','k')
if u.station_names
    text(x_sort,(x_sort.*0 )+ min(vp)-0.5,d.site(index),'rotation',u.station_names_angle,'interpreter','none','fontsize',8);
else % station names and title overlap. if station name shown, quanitity still shown in colorbar label
    title('Bahr Skew')
end

ylabel ('Log10 Period (s)')
axis([hp(1),hp(end), min(vp)-0.5 ,max(vp)])
caxis([0 0.4])
hcb = colorbar;
hcb.Label.String = 'Bahr Skew';
set(gca,'Layer','top')


annotation('textbox', [0 0.92 1 0.08], ...
'String', ['Bahr and Swift Skew Pseudo-Section | Profile azimuth = ',num2str(azimuth),char(176)], ...
'EdgeColor', 'none', ...
'HorizontalAlignment', 'center','FontSize',12,'FontWeight','bold')

