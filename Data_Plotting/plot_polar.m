function plot_polar(d,polar,is)
%%
% Function to plot polar diagram for a single site at 9 different periods
% pre-determined by the data
%
% Usage: plot_polar(d,polar,is)
%
% "d" is MT data structure
% "polar" is a structure of polar diagram variables
%           See calc_polar function
% "is" is site index
%
rot_ang = 0;

set_figure_size(1);

%Determine the 9 periods which will be plotted
[~, T_first]  = max(~isnan(d.Z(:,2,:)), [], 1); %First non-NaN period
T_first = min(T_first);

[~, T_last]  = max(flipud(~isnan(d.Z(:,2,:))), [], 1); %Last non-NaN period
T_last = d.nf-min(T_last);

n = 1;
for ifreq = floor(linspace(T_first,T_last,9)) %Loop over 9 periods

    subplot(3,3,n) %Will output a 3 x 3 subplot

    %Normalized polar diagram radius
    r = 1.1*max([squeeze(abs(polar.x(ifreq,2,is,:))); squeeze(abs(polar.y(ifreq,2,is,:)))]);

    % Plot whole ellipse
    plot(squeeze(polar.y(ifreq,2,is,:)),squeeze(polar.x(ifreq,2,is,:)),'k-') % Off-diagonal
    hold on
    plot(squeeze(polar.y(ifreq,1,is,:)),squeeze(polar.x(ifreq,1,is,:)),'k:') % Diagonal
    axis equal

    % Now plot axis showing coordinate frame
    xc1 = [0, r*cos(rot_ang*pi/180.)];      yc1 = [0, r*sin(rot_ang*pi/180.)] ;  
    xc2 = [0, r*cos(pi/2+rot_ang*pi/180.)]; yc2 = [0, r*sin(pi/2+rot_ang*pi/180.)];   
    plot(yc1,xc1,'r-'); plot(yc2,xc2,'b-');
    if ~isnan(r)
       axis([-r r -r r]);
       axis equal
       title(['T = ',num2str(d.T(ifreq)),' s'])
    else % if there is no data at this period
       title(['No Data for T = ',num2str(d.T(ifreq)),' s'])
    end

    n=n+1;
end

annotation('textbox', [0 0.92 1 0.08], ...
'String', ['Polar Diagram For Site ',d.site{is}], ...
'EdgeColor', 'none', ...
'HorizontalAlignment', 'center','FontSize',12,'FontWeight','bold')







