function plot_rho_pha(d,is)
%
% Function which plots resistivity and phase data for a single stations
%
% Usage: plot_rho_pha(d,is)
%
% "d" is an MT data structure
% "is" is the index of the station to plot (note: if "d" only includes one
% station then is must equal 1)
%
% If the data structure includes diagonal components, then the function
% will plot a 3 x 2 set of subplots which includes off-diagonal and
% diagonal plots as well as a map indicating which station is plotted. The
% bottom right subplot is left intentionally blank and is sometimes used to 
% plot other relevant things (e.g. rms)
%
% If there are no diagonal components, then the function will plot a 2 x 2
% set of subplots which includes only off-diagonal components and a map
% indicating which station is plotted. The bottom right subplot is left
% intentionally blank and is sometimes used to plot other relevant things
% (e.g. rms)

u = user_defaults;

if isnan(d.Zerr) %If all error is NaN then we are looking at predicted data and it will be plotted as a line
    linetype = {'-r','-b','-m','-g'};
else %If errors exist then we have observed data which will be plotted as points
    linetype = {'*r','*b','om','sg'};
end

% if isempty(find(ismember(d.responses,'ZXX'),1)) %If diagonals exist, set subplots accordingly
if ~strcmp('ZXX',d.responses)
    subplot(2,2,1);
else
    subplot(2,3,2);
end

%Plot off-diagonal apparent resistivities
logerrorbar(d.T,d.rho(:,2,is),d.rhoerr(:,2,is),linetype{1},'-r');hold on
logerrorbar(d.T,d.rho(:,3,is),d.rhoerr(:,3,is),linetype{2},'-b');
axis([u.Tlim(1) u.Tlim(2) u.rholim(1) u.rholim(2)])
ylabel('Apparent Resistivity (\Omega m)')
xlabel('Period (s)')
title(['Station ',num2str(is),' out of ',num2str(d.ns),': ',char(d.site(is))],'Interpreter','none');
grid on
manual_legend('XY',linetype{1},'YX',linetype{2});


if ~strcmp('ZXX',d.responses) %If diagonals exist, set subplots accordingly
    subplot(2,2,3);
else
    subplot(2,3,5);
end

%Plot off-diagonal phases
%Note: Phases are in e^{+iwt} convention so the yx phase is in the -180 to
%-90 quadrant. Adding 180 brings it into the correct quadrant
errorbar(d.T,d.pha(:,2,is),d.phaerr(:,2,is),linetype{1}); hold on
errorbar(d.T,d.pha(:,3,is)+180,d.phaerr(:,3,is),linetype{2});
set(gca, 'XScale', 'log')
axis([u.Tlim(1) u.Tlim(2) u.phalim(1) u.phalim(2)])
ylabel('Phase (deg)')
xlabel('Period (s)')
grid on

if any(strcmp('ZXX',d.responses))
% if isempty(find(ismember(d.responses,'ZXX'),1))~=1 %If diagonals exist, plot them
    
    %Plot diagonal apparent resistivities
    subplot(2,3,1); 
    logerrorbar(d.T,d.rho(:,1,is),d.rhoerr(:,1,is),linetype{3},'-m');hold on
    logerrorbar(d.T,d.rho(:,4,is),d.rhoerr(:,4,is),linetype{4},'-g');
    axis([u.Tlim(1) u.Tlim(2) u.rholimdiag(1) u.rholimdiag(2)])
    ylabel('Apparent Resistivity (\Omega m)')
    xlabel('Period (s)')
    grid on
    manual_legend('XX',linetype{3},'YY',linetype{4});

    %Plot diagonal phases
    %Note: Old way would add 180 deg to the yy phase. now just plot both
    %xx/yy phase in default [-180 180] since phase wraps are common
    subplot(2,3,4); 
    errorbar(d.T,d.pha(:,1,is),d.phaerr(:,1,is),linetype{3}); hold on
    errorbar(d.T,d.pha(:,4,is),d.phaerr(:,4,is),linetype{4});
    set(gca, 'XScale', 'log')
    axis([u.Tlim(1) u.Tlim(2) u.phalimdiag(1) u.phalimdiag(2)])
    ylabel('Phase (deg)')
    xlabel('Period (s)')
    grid on
end





end
