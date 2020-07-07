function plot_tipper(d,is,flag)
%
% Function which plots tipper for a single station
%
% Usage: plot_tipper(d,is,flag)
%
% "d" is an MT data structure
% "is" is the index of the station to plot (note: if "d" only includes one
% station then is must equal 1)
% "flag" is an OPTIONAL input 1 or 0 to set different subplot settings (0 is default)
%
% The function will plot a 2 x 2 set of subplots which includes real and
% imaginary Tx and Ty as well as a map indicating which station is plotted. The
% bottom right subplot is left intentionally blank and is sometimes used to 
% plot other relevant things (e.g. rms)

if ~exist('flag','var')
    flag = 0;
end

u = user_defaults;

if ~all(isnan(d.tip(:))) %If the entire tipper array is NaN, then no tipper data exists and do not plot anything
    

    if isnan(d.tiperr) %If all error is NaN then we are looking at predicted data and it will be plotted as a line
        linetype = {'-r','-b','-m','-g'};
    else %If errors exist then we have observed data which will be plotted as points
        linetype = {'or','sb','om','sg'};
    end

    %Plot real tipper components
    if flag
        subplot(2,3,3)
    else
        subplot(2,2,1)
    end
    errorbar(d.T,real(d.tip(:,1,is)),real(d.tiperr(:,1,is)),linetype{1});hold on
    errorbar(d.T,real(d.tip(:,2,is)),real(d.tiperr(:,2,is)),linetype{2});
    axis([u.Tlim(1) u.Tlim(2) u.tiplim(1) u.tiplim(2)])
    set(gca,'XScale','log')
    ylabel('Real Tipper')
    xlabel('Period (s)');
    title(['Station ',num2str(is),' out of ',num2str(d.ns),': ',char(d.site(is))],'Interpreter','none');
    h = zeros(2, 1);
    h(1) = plot(NaN,NaN,linetype{1});
    h(2) = plot(NaN,NaN,linetype{2});
    legend(h, 'Tx','Ty');

    %Plot imaginary tipper components
    if flag
        subplot(2,3,6)
    else
        subplot(2,2,3)
    end
    errorbar(d.T,imag(d.tip(:,1,is)),imag(d.tiperr(:,1,is)),linetype{1}); hold on
    errorbar(d.T,imag(d.tip(:,2,is)),imag(d.tiperr(:,2,is)),linetype{2});
    set(gca,'XScale','log')
    axis([u.Tlim(1) u.Tlim(2) u.tiplim(1) u.tiplim(2)])
    ylabel('Imaginary Tipper')
    xlabel('Period (s)')

else
    disp('No Tipper data in your observed data!')
end


end
