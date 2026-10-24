function plot_polar_ellipse_pseudo(d)
%
% Function which plots polar diagrams (peanuts) in pseudo-section format
% along a profile. Each polar diagram is plotted at the position of the
% station along the profile and at the height of its period.
%   - solid black line : off-diagonal impedance (Zxy)
%   - red line         : diagonal impedance (Zxx)
%
% Usage: plot_polar_ellipse_pseudo(d)
%
% "d" is an MT data structure
%
%%
close all
u = user_defaults;

% Size of the polar diagrams relative to the cell of the pseudo-section
% (1 = diagrams fill the cell). Optional: u.polar_pseudo_scale in user_defaults
if isfield(u,'polar_pseudo_scale')
    pscale = u.polar_pseudo_scale;
else
    pscale = 1;
end

%Pull a profile to plot the pseudo-section. "sidx" is stations indices
%included on profile and "azimuth" is the profile azimuth
[sidx,midpoint,azimuth] = get_pseudo_section_indices(d);

rad = 180./pi;

% Process cooordinates
lon_mean = midpoint(1);  lat_mean = midpoint(2);
x = cos(lat_mean/rad)*111*(d.loc(sidx,2)-lon_mean);y = 111*(d.loc(sidx,1)-lat_mean); % km
N = length(sidx); %Number of stations included in profile

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

% Polar diagram coordinates for all stations / periods
% polar.x, polar.y have size (nf x 4 x ns x nrot+1)
[polar] = calc_polar(d);

%%
%Plot pseudo-section
set_figure_size(1);

logT     = log10(d.T(:));
rangeT   = max(logT)-min(logT);
dlogT    = u.nskip*rangeT/max(d.nf-1,1);   % spacing (log10 period) between plotted periods

% Horizontal position of each station and aspect ratio of the section
if u.profile_distance
    xpos = x_sort(:);                       % true distance along profile (km)
    if rangeT == 0 || max(xpos)-min(xpos) == 0 % 1 station or 1 period
        aspect_ratio = 1;
    else
        aspect_ratio = (max(xpos)-min(xpos))/rangeT;
    end
else
    xpos = (1:N)';                          % equal spacing between stations
    if dlogT == 0
        aspect_ratio = 1;
    else
        aspect_ratio = 1/dlogT;             % one period step looks as large as one station step
    end
end

% Size of one cell of the pseudo-section in x-axis units
if N > 1
    sx = median(diff(xpos));                % station spacing
else
    sx = aspect_ratio*dlogT;
end
sy = aspect_ratio*dlogT;                    % period spacing (in x-axis units)
if sy == 0; sy = sx; end
Rmax = 0.45*pscale*min(sx,sy);              % maximum radius of a polar diagram

for ifreq = 1:u.nskip:d.nf
    for is = 1:N %Loop of stations included on profile

        % Normalized polar diagram radius (based on the off-diagonal term)
        r = 1.1*max([squeeze(abs(polar.x(ifreq,2,index(is),:))); squeeze(abs(polar.y(ifreq,2,index(is),:)))]);
        if isnan(r) || r == 0
            continue % no data at this period / station
        end
        scale = Rmax/r;

        % X location on profile, and Y location in period.
        % Note the minus sign for y - this is to negate the flip from axis ij later on
        xoff = xpos(is) + scale*squeeze(polar.y(ifreq,:,index(is),:));  % 4 x nrot+1
        yoff = logT(ifreq) - scale*squeeze(polar.x(ifreq,:,index(is),:))./aspect_ratio;

        plot(xoff(2,:),yoff(2,:),'k-'); hold on % Off-diagonal
        plot(xoff(1,:),yoff(1,:),'r-');         % Diagonal

    end
end

daspect(gca,[aspect_ratio 1 1]);
axis ij
ylabel('Log(Period (s))')

if u.profile_distance

    set(gca,'fontsize',20)
    xlabel('Distance Along Profile (km)');

    %Plot station locations on profile
    dx = 0.05*(x_sort(end) - x_sort(1));
    plot(x_sort,min(logT)-0.25,'kv','MarkerFaceColor','k')

    if dx == 0
        axis([x_sort(1)-1,x_sort(N)+1, min(logT)-0.5 ,max(logT)+0.5*dlogT])
    else
        axis([x_sort(1)-dx,x_sort(N)+dx, min(logT)-0.5 ,max(logT)+0.5*dlogT])
    end

    if u.station_names
        text(x_sort,(x_sort.*0 )+ min(logT)-0.5,d.site(index),'rotation',u.station_names_angle,'interpreter','none');
    else
        title('Polar Diagram Pseudo Section')
    end

else

    xlabel('Station Name')

    h = gca;
    h.TickLabelInterpreter = 'none';

    xticks(xpos);
    set(h,'XTickLabel',d.site(index),'FontSize',12)

    axis([0.5, N+0.5, min(logT)-0.5 ,max(logT)+0.5*dlogT])

end

set(gca,'Layer','top');

print_figure('polar',['polar_pseudo_ellipse_',num2str(lat_mean),'_lat_',num2str(lon_mean),'_lon']); %Save figure

end