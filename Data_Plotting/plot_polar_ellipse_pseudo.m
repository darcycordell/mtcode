function plot_polar_ellipse_pseudo(d)
%
% Function which plots phase tensors ellipes in pseudo-section format and
% colors the ellipses with beta skew angle
%
% Usage: plot_phase_tensor_ellipse_pseudo(d)
%
% "d" is an MT data structure
%
%%
close all
u = user_defaults;

%Pull a profile to plot the pseudo-section. "sidx" is stations indices
%included on profile and "rot_ang" is the direction of the normal to the
%profile direction
[sidx,midpoint,azimuth] = get_pseudo_section_indices(d);

rad = 180./pi;      

% Process cooordinates
lon_mean = midpoint(1);  lat_mean = midpoint(2); % km
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

%Initialize arrays
[polar] = calc_polar(d);


%%
%Plot pseudo-section
set_figure_size(1);
%Set vertical scaling for periods and phase tensor ellipses
aspect_ratio = (max(x_sort)-min(x_sort))/(max(log10(d.T))-min(log10(d.T)));
if aspect_ratio == 0 % 1 station
    aspect_ratio = 1;
end

for ifreq = 1:u.nskip:d.nf
    for is = 1:N %Loop of stations included on profile

        r = 1.1*max([squeeze(abs(polar.x(ifreq,2,index(is),:))); squeeze(abs(polar.y(ifreq,2,index(is),:)))]);

    % Plot whole ellipse
    plot(squeeze(polar.y(ifreq,2,is,:)),squeeze(polar.x(ifreq,2,is,:)),'k-') % Off-diagonal
    hold on
    plot(squeeze(polar.y(ifreq,1,is,:)),squeeze(polar.x(ifreq,1,is,:)),'k:') % Diagonal
    axis equal

        if u.profile_distance
            scale = u.pt_pseudo_scale*((max(x_sort)-min(x_sort))/length(x_sort))/r*(((max(log10(d.T)))-min(log10(d.T)))/(d.nf));
            if scale == 0 % 1 station
                scale = u.pt_pseudo_scale/r*(((max(log10(d.f)))-min(log10(d.f)))/(d.nf));
            end
            xp = x_sort(is)+scale*polar.y(ifreq,index(is)); %X location on profile 
            %xp = 2*is+scale*pty(:,ifreq,index(is)); %Option to plot site-by-site rather as distance

        else
            xref = (1/(((max(log10(d.T)))-min(log10(d.T)))/(d.nf)));
            scale = u.pt_pseudo_scale*((max(x_sort)-min(x_sort))/length(x_sort))/r*(1/xref);

            xp = 1.2*xref*index(is) + scale*pty(:,ifreq,index(is));

            xpref(is) = 1.2*xref*index(is);

        end
            
            
        % Note the minus sign below - this is to negate the flip from axis ij later on
        yp = log10(d.T(ifreq))-scale*ptx(:,ifreq,index(is))./aspect_ratio;  %Y location in period
        if strcmp(u.phase_tensor_ellipse_fill,'phimin')
            fill(xp,yp,abs(phimin(ifreq,index(is)))); hold on %Plot filled ellipses with phi min
        elseif strcmp(u.phase_tensor_ellipse_fill,'beta')
            fill(xp,yp,abs(beta(ifreq,index(is)))); hold on %Plot filled ellipses with beta skew angle
        else
            disp('Unrecognized input for u.phase_tensor_ellipse_fill. Ellipses are filled with beta skew angle values. Check your user_defaults')
            fill(xp,yp,abs(beta(ifreq,index(is)))); hold on %Plot filled ellipses with beta skew angle
        end

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
    plot(x_sort,log10(d.T(1))-0.25,'kv','MarkerFaceColor','k')

    if dx ==0
        axis([x_sort(1)-1,x_sort(N)+1, min(log10(d.T))-0.5 ,max(log10(d.T))])
    else
        axis([x_sort(1)-dx,x_sort(N)+dx, min(log10(d.T))-0.5 ,max(log10(d.T))])
    end

    if u.station_names
        text(x_sort,(x_sort.*0 )+ min(log10(d.T))-0.5,d.site(index),'rotation',u.station_names_angle,'interpreter','none');
    else
        title('Phase Tensor Pseudo Section')
    end
    
else
    
    xlabel('Station Name')
    
    h = gca;
    h.TickLabelInterpreter = 'none';
    
    xticks([xpref]);
    set(h,'XTickLabel',d.site(index),'FontSize',12)
    
    axis([min(xpref)-0.75*min(xpref),max(xpref)+0.75*min(xpref), min(log10(d.T))-0.5 ,max(log10(d.T))])


end

set(gca,'Layer','top');


if strcmp(u.phase_tensor_ellipse_fill,'phimin')
    caxis(u.phase_tensor_phimin_colim);
    hcb = colorbar;
    hcb.Label.String = '\Phi_{min} (degrees)';
else
    caxis(u.phase_tensor_beta_colim);
    hcb = colorbar;
    hcb.Label.String = 'Beta Skew Angle (degrees)';
end



print_figure(['phase_tensor_',d.niter],['pseudo_ellipse_',num2str(lat_mean),'_lat_',num2str(lon_mean),'_lon']); %Save figure

