function plot_phase_tensor_pseudo(d)
%
% Function which plots pseudo-section of phase tensor parameters such as
% phimin, phimax, alpha, beta, eccentricity
%
% Usage: plot_phase_tensor_pseudo(d)
%
% "d" is an MT data structure
%
%%
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

% Put stations in order of monotonically increasing x_rot
[x_sort,index] = sort(x_rot,'ascend');
index = sidx(index);

%Initialize arrays
beta = nan(d.nf,d.ns); alpha = nan(size(beta)); e = nan(size(beta));
strike = nan(size(beta)); phimin = nan(size(beta)); phimax = nan(size(beta));

for is = 1:d.ns  % Loop over stations

    [p] = calc_phase_tensor(d.Z(:,:,is)); %Calculate phase tensor parameters

    %Put all the phase tensor parameters into matrices of size nf x ns
    beta(:,is) = p.beta;
    alpha(:,is) = p.alpha;
    e(:,is) = p.e;
    strike(:,is) = p.strike;
    phimin(:,is) = p.phimin;
    phimax(:,is) = p.phimax;

end

%========================================================================
%Plot alpha angle pseudo section
set_figure_size(1);
subplot(3,2,1)
% pcolor(x_sort,log10(d.T),alpha(:,index));
[hp,vp,C] = fix_pcolor(x_sort,log10(d.T),alpha(:,index));
pcolor(hp,vp,C)
hold on; axis ij
shading 'flat'
plot(x_sort,log10(d.T(1))-0.25,'kv','MarkerFaceColor','k')
if u.station_names
    text(x_sort,(x_sort.*0 )+ min(vp)-0.5,d.site(index),'rotation',u.station_names_angle,'interpreter','none','fontsize',8);
else
    title('Alpha ')
end

ylabel ('Log10 Period (s)')
% axis([x_sort(1),x_sort(N), min(log10(d.T))-0.5 ,max(log10(d.T))])
axis([hp(1),hp(end), min(vp)-0.5 ,max(vp)])
hcb = colorbar;
hcb.Label.String = 'Alpha Angle (degrees)';
set(gca,'Layer','top')

%========================================================================
%Plot beta skew angle
subplot(3,2,2)
[hp,vp,C] = fix_pcolor(x_sort,log10(d.T),abs(beta(:,index)));
pcolor(hp,vp,C)
hold on; axis ij
shading 'flat'
plot(x_sort,log10(d.T(1))-0.25,'kv','MarkerFaceColor','k')
if u.station_names
    text(x_sort,(x_sort.*0 )+ min(vp)-0.5,d.site(index),'rotation',u.station_names_angle,'interpreter','none','fontsize',8);
else
    title('Beta ')
end

ylabel ('Log10 Period (s)')
axis([hp(1),hp(end), min(vp)-0.5 ,max(vp)])
caxis(u.phase_tensor_beta_colim)
hcb = colorbar;
hcb.Label.String = 'Beta Skew Angle (degrees)';
set(gca,'Layer','top')

%========================================================================
%Plot minimum phase tensor axis length
subplot(3,2,3)
[hp,vp,C] = fix_pcolor(x_sort,log10(d.T),phimin(:,index));
pcolor(hp,vp,C)
hold on; axis ij
shading 'flat'
plot(x_sort,log10(d.T(1))-0.25,'kv','MarkerFaceColor','k')
if u.station_names
    text(x_sort,(x_sort.*0 )+ min(vp)-0.5,d.site(index),'rotation',u.station_names_angle,'interpreter','none','fontsize',8);
else
    title('Minimum Phase Tensor Ellipse Axis')
end

ylabel ('Log10 Period (s)')
axis([hp(1),hp(end), min(vp)-0.5 ,max(vp)])
hcb = colorbar;
hcb.Label.String = 'Length of Ellipse Axis';
set(gca,'Layer','top')

%========================================================================
%Plot maximum phase tensor axis length
subplot(3,2,4)
[hp,vp,C] = fix_pcolor(x_sort,log10(d.T),phimax(:,index));
pcolor(hp,vp,C)
hold on; axis ij
shading 'flat'
plot(x_sort,log10(d.T(1))-0.25,'kv','MarkerFaceColor','k')
if u.station_names
    text(x_sort,(x_sort.*0 )+ min(vp)-0.5,d.site(index),'rotation',u.station_names_angle,'interpreter','none','fontsize',8);
else
    title('Maximum Phase Tensor Ellipse Axis')
end

ylabel ('Log10 Period (s)')
axis([hp(1),hp(end), min(vp)-0.5 ,max(vp)])
hcb = colorbar;
hcb.Label.String = 'Length of Ellipse Axis';
set(gca,'Layer','top')

%========================================================================
%Plot phase tensor eccentricity
subplot(3,2,5)
[hp,vp,C] = fix_pcolor(x_sort,log10(d.T),e(:,index));
pcolor(hp,vp,C)
hold on; axis ij
shading 'flat'
plot(x_sort,log10(d.T(1))-0.25,'kv','MarkerFaceColor','k')
if u.station_names
    text(x_sort,(x_sort.*0 )+ min(vp)-0.5,d.site(index),'rotation',u.station_names_angle,'interpreter','none','fontsize',8);
else
    title('Phase Tensor Eccentricity')
end

ylabel ('Log10 Period (s)')
axis([hp(1),hp(end), min(vp)-0.5 ,max(vp)])
hcb = colorbar;
hcb.Label.String = 'Eccentricity';
set(gca,'Layer','top')

%========================================================================
%Plot phase tensor geoelectric strike (alpha - beta)
subplot(3,2,6)
[hp,vp,C] = fix_pcolor(x_sort,log10(d.T),strike(:,index));
pcolor(hp,vp,C)
hold on; axis ij
shading 'flat'
plot(x_sort,log10(d.T(1))-0.25,'kv','MarkerFaceColor','k')
if u.station_names
    text(x_sort,(x_sort.*0 )+ min(vp)-0.5,d.site(index),'rotation',u.station_names_angle,'interpreter','none','fontsize',8);
else % station names and title overlap. if station name shown, quanitity still shown in colorbar label
    title('Geoelectric Strike')
end

ylabel ('Log10 Period (s)')
axis([hp(1),hp(end), min(vp)-0.5 ,max(vp)])
hcb = colorbar;
hcb.Label.String = 'Strike Angle';
set(gca,'Layer','top')

annotation('textbox', [0 0.92 1 0.08], ...
'String', ['Phase Tensor Pseudo-Section | Profile azimuth = ',num2str(azimuth),char(176)], ...
'EdgeColor', 'none', ...
'HorizontalAlignment', 'center','FontSize',12,'FontWeight','bold')

print_figure(['phase_tensor_',d.niter],['pseudo_',num2str(lat_mean),'_lat_',num2str(lon_mean),'_lon']); %Save figure

