function plot_phase_tensor_beta_strike(d)
%%
% Function which plots a summary of phase tensor parameters (beta angle,
% geoelectric strike, eccentricity) as a function of period for all
% stations. The mean of the parameters is computed for all stations.
%
% Also plots rose diagram histograms of strike and histograms of beta skew
% and eccentricity at each period.
%
% Usage: plot_phase_tensor_beta_strike(d)
%
% "d" is the input data structure.

u = user_defaults;

beta = nan(d.nf,d.ns);
alpha = nan(size(beta));
e = nan(size(beta));
strike = nan(size(beta));
ratio = nan(size(beta));
for is = 1:d.ns  % Loop over stations

    [p] = calc_phase_tensor(d.Z(:,:,is)); %Calculate phase tensor info
        %The "p" structure contains all info for phase tensor

    %Put all the outputs into matrices of size nf x ns
    beta(:,is) = p.beta;
    alpha(:,is) = p.alpha;
    e(:,is) = p.e;
    strike(:,is) = p.strike;
    ratio(:,is) = abs(p.phimax./p.phimin);

end


for ifreq = 1:u.nskip:d.nf  % Loop over frequencies

    set_figure_size(1);
    
    subplot(2,2,1)
    strike_vec = [strike(ifreq,:) strike(ifreq,:)+90 strike(ifreq,:)+180 strike(ifreq,:)+270];
    rose_geog(strike_vec,24,u.rose_histogram,'b')
    
    [avg_strike, std] = strike_stats(strike_vec,4);

    title(['Avg \alpha - \beta = ', num2str(avg_strike),'^{\circ}; \sigma = ',num2str(std),'^{\circ}'])


    subplot(2,2,2)
    strike_vec = [alpha(ifreq,:) alpha(ifreq,:)+90 alpha(ifreq,:)+180 alpha(ifreq,:)+270];
    rose_geog(strike_vec,24,u.rose_histogram,'b')

    if d.T(ifreq) > 1
      title(['\alpha : T = ',num2str(d.T(ifreq)),' s'])
    else
      title(['\alpha : f = ',num2str(1/d.T(ifreq)),' Hz'])
    end


    subplot(2,2,3)
    xint=0.1;
    edges=1:xint:2;
    [rn,~]=histc(ratio(ifreq,:),edges);
    bar(edges,rn,'histc');
    title(['median phi ratio = ',num2str( median(ratio(ifreq,~isnan(e(ifreq,:)))))])
    xlabel('\phi_{max} / \phi_{min}'); ylabel('Number of Stations');
    axis([1 2 0 u.rose_histogram]);

    subplot(2,2,4)
    xint=0.5;
    edges=0:xint:5;
    [betan,~]=histc(abs(beta(ifreq,:)),edges);
    bar(edges,betan,'histc');
    xlabel('\beta skew'); ylabel('Number of Stations');
    axis([0 5 0 u.rose_histogram]);
    title(['median beta skew = ',num2str(median(abs(beta(ifreq,~isnan(e(ifreq,:))))))])
    
    pause(0.1)
    
    print_figure(['phase_tensor_',d.niter],['phase_tensor_rose_',num2str(ifreq,'%03.0f'),'_',num2str(d.T(ifreq)),'_s']); %Save figure


end



 
%Plot the mean alpha angle
set_figure_size(1);
subplot(4,1,1)
semilogx(d.T,nanmean(alpha,2),'LineWidth',1.5)
ylabel('Alpha (deg)'); grid on
title('Average Phase Tensor Parameters vs. Period')

%Plot the mean beta skew angle
subplot(4,1,2)
semilogx(d.T,nanmean(abs(beta),2),'LineWidth',1.5);
ylabel('Beta Skew (deg)'); grid on

%Plot the mean geoelectric strike angle (alpha - beta)
subplot(4,1,3)
semilogx(d.T,nanmean(strike,2),'LineWidth',1.5);
ylabel('Strike \alpha - \beta (deg)'); grid on

%Plot the mean eccentricity
subplot(4,1,4)
semilogx(d.T,nanmean(e,2),'LineWidth',1.5);
xlabel('Period (s)'); grid on
ylabel('Eccentricity')

print_figure(['phase_tensor_',d.niter],'phase_tensor_summary'); %Save figure




end



