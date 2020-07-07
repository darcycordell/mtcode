function plot_impedance(d,is)
% Function which plots impedance data for a single station
%
% Usage: plot_impedance(d,is)
%
% "d" is MT data structure
% "is" is the site index to plot
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
%

u = user_defaults;

if isnan(d.Zerr) %If all error is NaN then we are looking at predicted data and it will be plotted as a line
    linetype = {'-r','-b','-m','-g'};
else %If errors exist then we have observed data which will be plotted as points
    linetype = {'or','sb','om','sg'};
end

if isempty(find(ismember(d.responses,'ZXX'),1)) %If diagonals exist, set subplots accordingly
    subplot(2,2,1)
else
    subplot(2,3,2)
end

%Plot real off-diagonal impedances
logerrorbar(d.T,abs(real(d.Z(:,2,is))),real(d.Zerr(:,2,is)),linetype{1},'-r');hold on
logerrorbar(d.T,abs(real(d.Z(:,3,is))),real(d.Zerr(:,3,is)),linetype{2},'-b');
axis([u.Tlim(1) u.Tlim(2) u.Zlim(1) u.Zlim(2)])
ylabel('Real Off-Diagonal Impedance')
xlabel('Period (s)');
title(['Station ',num2str(is),' out of ',num2str(d.ns),': ',char(d.site(is))],'Interpreter','none');
grid on
manual_legend('XY',linetype{1},'YX',linetype{2});


if isempty(find(ismember(d.responses,'ZXX'),1)) %If diagonals exist, set subplots accordingly
    subplot(2,2,3)
else
    subplot(2,3,5)
end
%Plot imaginary off-diagonal impedances
logerrorbar(d.T,abs(imag(d.Z(:,2,is))),imag(d.Zerr(:,2,is)),linetype{1},'-r'); hold on
logerrorbar(d.T,abs(imag(d.Z(:,3,is))),imag(d.Zerr(:,3,is)),linetype{2},'-b');
axis([u.Tlim(1) u.Tlim(2) u.Zlim(1) u.Zlim(2)])
ylabel('Imaginary Off-Diagonal Impedance')
xlabel('Period (s)')
grid on


if isempty(find(ismember(d.responses,'ZXX'),1))~=1 %If diagonals exist, plot them
    subplot(2,3,1)
    %Plot real diagonal impedances
    logerrorbar(d.T,abs(real(d.Z(:,1,is))),real(d.Zerr(:,1,is)),linetype{3},'-m');hold on
    logerrorbar(d.T,abs(real(d.Z(:,4,is))),real(d.Zerr(:,4,is)),linetype{4},'-g');
    axis([u.Tlim(1) u.Tlim(2) u.Zlim(1) u.Zlim(2)])
    ylabel('Real Diagonal Impedance')
    xlabel('Period (s)')
    grid on
    manual_legend('XX',linetype{3},'YY',linetype{4});

    %Plot imaginary diagonal impedances
    subplot(2,3,4)
    logerrorbar(d.T,abs(imag(d.Z(:,1,is))),imag(d.Zerr(:,1,is)),linetype{3},'-m'); hold on
    logerrorbar(d.T,abs(imag(d.Z(:,4,is))),imag(d.Zerr(:,4,is)),linetype{4},'-g');
    axis([u.Tlim(1) u.Tlim(2) u.Zlim(1) u.Zlim(2)])
    ylabel('Imaginary Off-Diagonal Impedance')
    xlabel('Period (s)')
    grid on
    
end




end
