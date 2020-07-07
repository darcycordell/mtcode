function plot_misfit_rms_by_period(dobs,dpred,s)
% Function which plots rms misfit between two MT data sets as a function of
% period for all sites, overall and by impedance component.
%
% Usage: 
% plot_misfit_rms_by_period(dobs,dpred)
% plot_misfit_rms_by_period(dobs,dpred,s)
%
% "dobs" is observed MT data. Note: the errors used to normalize the
% residuals are taken from dobs.
% "dpred" is the predicted MT data. The errors in dpred are ignored.
% "s" is an OPTIONAL input that is the stats structure output from
% detailed_statistics or detailed_staistics_2D. If s is not an input, the
% misfit stats are computed from impedances by default.
%
%%
u = user_defaults;

%Calculate detailed statistics for all data and also rms by frequency, rms
%by site by frequency, and rms by frequency by component
if ~exist('s','var')
    disp('No stats structure input. Calculating rms misfit from impedances')
    [s] = detailed_statistics(dobs,dpred);
end
rms_freq = s.rms_freq;
rms_sf = s.rms_sf;
rms_fr = s.rms_fr;

close all

for i = 1:dobs.ns %Plot rms by period for all sites individually
    h(1) = semilogx(dobs.T,rms_sf(i,:),'-b'); hold on
    h(1).Color(4) = 0.25;
end

%Plot rms by period overall
h(2) = semilogx(dobs.T,rms_freq,'-r','LineWidth',3); hold on

%Plot overall rms as a straight line on the plot
h(3) = semilogx([min(dobs.T) max(dobs.T)],[s.rms s.rms],'--k');
xlabel('Period (s)')
ylabel('Root Mean Square Misfit')
axis([u.Tlim(1) u.Tlim(2) 0 max(1.1*rms_sf(:))])
legend(h,'Individual Stations','By Period','Overall RMS')
print_figure(['misfit_stats_',dpred.niter],'rms_by_period'); %Save figure


figure(2) %Plot rms by period by component
if isempty(find(ismember(dobs.responses,'TX'),1)) %If tipper exists, set up subplots accordingly
    subplot(2,1,1)
else
    subplot(3,1,1)
end
%Plot rms for diagonal components by period
semilogx(dobs.T,rms_fr(1,:),'o-m'); hold on
semilogx(dobs.T,rms_fr(4,:),'.-g');
xlabel('Period (s)')
ylabel('Root Mean Square Misfit')
legend('XX','YY')
axis([u.Tlim(1) u.Tlim(2) 0 nanmax([rms_fr(:); 5])])
title('Root Mean Square: Diagonal Impedances')

if isempty(find(ismember(dobs.responses,'TX'),1)) %If tipper exists, set up subplots accordingly
    subplot(2,1,2)
else
    subplot(3,1,2)
end
%Plot rms for off-diagonal components by period
semilogx(dobs.T,rms_fr(2,:),'o-r'); hold on
semilogx(dobs.T,rms_fr(3,:),'.-b');
xlabel('Period (s)')
ylabel('Root Mean Square Misfit')
legend('XY','YX')
axis([u.Tlim(1) u.Tlim(2) 0 nanmax([rms_fr(:); 5])])
title('Root Mean Square: Off-Diagonal Impedances')


if isempty(find(ismember(dobs.responses,'TX'),1))~=1 %If tipper exists, plot them
    subplot(3,1,3)
    %Plot rms for tipper components by period
    semilogx(dobs.T,rms_fr(5,:),'o-k');hold on
    semilogx(dobs.T,rms_fr(6,:),'.-c');
    xlabel('Period (s)')
    ylabel('Root Mean Square Misfit')
    legend('Tx','Ty')
    axis([u.Tlim(1) u.Tlim(2) 0 max(rms_fr(:))])
    title('Root Mean Square: Tipper')
end

print_figure(['misfit_stats_',dpred.niter],'rms_by_period_by_component'); %Save figure
