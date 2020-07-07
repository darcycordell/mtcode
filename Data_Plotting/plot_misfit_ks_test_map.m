function [ks,x_all,y_all] = plot_misfit_ks_test_map(dtrue,dpred1,dpred2,alpha,flag)
%
% Function which calculates the Kolmogorov-Smirnov statistical test (KS
% test) for two distributions of residuals. The KS test is performed for
% the entire dataset as well as station-by-station. The station-by-station
% results are plotted in map view.
%
% The KS test tests whether two distributions are statistically equivalent
% (e.g. are the distributions the same with a 95% confidence?)
%
% The KS test returns a p-value. If the p-value is less than the confidence
% level (e.g. alpha = 0.05), then there is a "statistically significant"
% difference between the two distributions.
%
% For MT, we might have two distributions of residuals. One set of
% residuals came from the comparison between an inversion and true data,
% and the other set of residuals from the comparison between a different
% inversion and true data. We always compare the residuals to the true
% data. Then we look at if the distribution of residuals is different.
%
% A menu will ask if you would like to use impedance data, tipper data, or
% both impedance and tipper data.
%
% Note that the KS test works even if the length of the two vectors is
% different. So you can compare different data sets with different numbers
% of stations, periods, etc. However, of course you cannot make cross-plots of these
% datasets.
%
% Usage: [ks] = plot_misfit_ks_test_map(dtrue,dpred1,dpred2,alpha,flag)
%
% Inputs: dtrue is a "d" data structure of "true" or "real" MT data
%       dpred1 is a data structure of predicted data (can be forward or inversion response)
%       dpred2 is another data structure of predicted data
%       alpha is the confidence level (default is 0.05)
%       flag is an option argument (1 or 0) to plot or not plot
%
% Outputs: ks is a structure which contains a p-value for all the residuals
% and a true-false value showing whether the distributions are different
% (h=1) or the same (h=0). The ks structure also incldues a vector p-values
% for each MT site.
%
% 
% 
%%
u = user_defaults;
[L] = load_geoboundary_file_list;

close all

if ~exist('flag','var')
    flag = 1;
end

%Get residuals for each predicted dataset compared to the true dataset
[s1] = detailed_statistics(dtrue,dpred1);
[s2] = detailed_statistics(dtrue,dpred2);

response_menu = menu('Choose data type for K-S test','Impedance + Tipper (if existing)','Impedance only','Tipper only');

    %Make a vector of residuals for each dataset
if response_menu == 1 % impedance and tipper
    x_all = mv(real(s1.residuals.imp),imag(s1.residuals.imp),real(s1.residuals.tip),imag(s1.residuals.tip)); % residuals between dtrue and dpred1
    y_all = mv(real(s2.residuals.imp),imag(s2.residuals.imp),real(s2.residuals.tip),imag(s2.residuals.tip)); % residuals between dtrue and dpred2

elseif response_menu == 2 % impedance only
    x_all = mv(real(s1.residuals.imp),imag(s1.residuals.imp)); % residuals between dtrue and dpred1
    y_all = mv(real(s2.residuals.imp),imag(s2.residuals.imp)); % residuals between dtrue and dpred2

elseif response_menu == 3 % tipper only
    x_all = mv(real(s1.residuals.tip),imag(s1.residuals.tip)); % residuals between dtrue and dpred1
    y_all = mv(real(s2.residuals.tip),imag(s2.residuals.tip)); % residuals between dtrue and dpred2
    
else
    return
end

x_all(isnan(x_all))=[];
y_all(isnan(y_all))=[];

[h_all,p_all,~] = kstest2(x_all,y_all,'Alpha',alpha); %KS test (built in MATLAB function)

% if flag && length(x_all)==length(y_all) %Make cross-plot only if the length of the two vectors is the same
    %Plot all residuals. If the two distributions are the same, it should plot
    %as a straight line.
%     figure(1)
%     plot(x_all,y_all,'.k'); hold on
%     plot([-1000 1000],[-1000 1000],'--r');
%     axis([min(x_all) max(x_all) min(y_all) max(y_all)])
%     xlabel('Dataset #1 Normalized Residuals')
%     ylabel('Dataset #2 Normalized Residuals')
%     title('All Residuals')
%     plot_misfit_cross_plot(dtrue,dpred1,dpred2,1);
%     print_figure(['statistics_compare'],['ks_test_cross_plot']); %Save figure
% end

%Compute the KS test for each station only if the dataset contains the same
%number of stations
if dpred1.ns == dpred2.ns
    for i = 1:dtrue.ns
        
        if response_menu == 1 % impedance and tipper
            x = mv(real(s1.residuals.imp(:,:,i)),imag(s1.residuals.imp(:,:,i)),real(s1.residuals.tip(:,:,i)),imag(s1.residuals.tip(:,:,i)));
            y = mv(real(s2.residuals.imp(:,:,i)),imag(s2.residuals.imp(:,:,i)),real(s2.residuals.tip(:,:,i)),imag(s2.residuals.tip(:,:,i)));

        elseif response_menu == 2 % impedance only
            x = mv(real(s1.residuals.imp(:,:,i)),imag(s1.residuals.imp(:,:,i)));
            y = mv(real(s2.residuals.imp(:,:,i)),imag(s2.residuals.imp(:,:,i)));
        
        elseif response_menu == 3 % tipper only
            x = mv(real(s1.residuals.tip(:,:,i)),imag(s1.residuals.tip(:,:,i)));
            y = mv(real(s2.residuals.tip(:,:,i)),imag(s2.residuals.tip(:,:,i)));
        
        end

        x(isnan(x)) = [];
        y(isnan(y)) = [];

        %[p,h(i)] = ranksum(x,y,'alpha',0.05); %Other statistical test, I think
        %KS is better
        if isempty(x) || isempty(y)
            h(i) = false;
            p(i) = NaN;
        else
            [h(i),p(i)] = kstest2(x,y,'Alpha',alpha);
        end

    %     figure
    %     plot(x,y,'.k'); hold on
    %     plot([-100 100],[-100 100],'--r');
    %     axis([min(x) max(x) min(y) max(y)])
    %     axis equal
    %     title(num2str(h))

    end

    if flag
        % %Plot KS test p-values for each site in map view
        set_figure_size(2);
        for i = 1:dtrue.ns %Plot each site rms as a colored circle at the site location in map view

            th = 0:pi/20:2*pi;
            r = abs(dtrue.lim(4)-dtrue.lim(3))/50; %This is the radius of the circle. It is 1/50 of the survey area
            xp = r * cos(th) + dtrue.loc(i,2);
            yp = r * sin(th)*0.85+ dtrue.loc(i,1);

            m_fill(xp,yp,p(i)); hold on;
            shading flat;

        end

        plot_geoboundaries(L);

        if h_all == 1
            pass = 'Pass';
        else
            pass = 'Fail';
        end

        caxis([0 1]);
        colormap(u.cmap)
        hcb = colorbar;
        hcb.Label.String = 'P Value';
        hcb.FontSize = 12;
        m_grid('box','fancy','tickdir','in','xlabeldir','end'); hold on
        title(['KS Test ',pass,'. Overall P-value = ',num2str(p_all),'. ',num2str(length(find(h==1))),' Sites with "Significant" Change']);


%         print_figure('statistics_compare',['ks_test_pval_',num2str(p_all)]); %Save figure

    end

end

ks.p = p;
ks.p_all = p_all;
ks.h = h;
ks.h_all = h_all;

