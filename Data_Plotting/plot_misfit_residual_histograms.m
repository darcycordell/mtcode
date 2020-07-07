function plot_misfit_residual_histograms(varargin)
% Function which plots residuals between two data set structures on a
% histogram. Plots overall residuals as well as residuals for each component
%
% Usage: plot_misfit_residual_histograms(varargin)
%
% Can use this function with 1 or 2 inputs:
%
% 1 input :
% plot_misfit_residual_histograms(s)
% "s" is a stats structure output from detailed_statistics or
% detailed_statistics_2D
%
% 2 inputs :
% plot_misfit_residual_histograms(dobs,dpred)
% "dobs" is observed MT data structure
% "dpred" is predicted MT data structure
%
% If a component does not exist, the plot is blank.

%Calculate detailed statistics for all data and also rms by component
if nargin == 2 % dobs and dpred are input
    disp('Observed and predicted data are inputs. Calculating impedance residuals')
    [s] = detailed_statistics(varargin{1},varargin{2});
elseif nargin == 1 % only s structure input
    s = varargin{1};
else
    error('Unrecognized number of inputs')
end

close all
set_figure_size(1);

nbins = 500;
    
rms_comp = s.rms_comp;
residuals = s.residuals;

if isfield(s.residuals,'rho')
    
    disp('Plotting residuals of apparent resistivity and phase')
    
    %Plot overall residual histograms
    subplot(3,3,[1 4 7])
    histfit([real(residuals.tip(:)); imag(residuals.tip(:)); residuals.rho(:); residuals.pha(:)],nbins); hold on
    xlabel('Normalized Residual'); ylabel('Number of Data Points')
    title(['All Datapoints: RMS = ',num2str(s.rms)])

    subplot(3,3,2) %XX component
    histfit(mv((residuals.rho(:,1,:)),(residuals.pha(:,1,:))),nbins)
    xlabel('Normalized Residual'); ylabel('Number of Data Points')
    title(['Rho and Phase XX: RMS = ',num2str(rms_comp(1))])  % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<check components in 2D

    subplot(3,3,3) %%XY component
    histfit(mv((residuals.rho(:,2,:)),(residuals.pha(:,2,:))),nbins)
    xlabel('Normalized Residual'); ylabel('Number of Data Points')
    title(['Rho and Phase XY: RMS = ',num2str(rms_comp(2))])

    subplot(3,3,5) %YX component
    histfit(mv((residuals.rho(:,3,:)),(residuals.pha(:,3,:))),nbins)
    xlabel('Normalized Residual'); ylabel('Number of Data Points')
    title(['Rho and Phase YX: RMS = ',num2str(rms_comp(3))])

    subplot(3,3,6) %%YY component
    histfit(mv(real(residuals.rho(:,4,:)),imag(residuals.pha(:,4,:))),nbins)
    xlabel('Normalized Residual'); ylabel('Number of Data Points')
    title(['Rho and Phase YY: RMS = ',num2str(rms_comp(4))])

    subplot(3,3,8) %TX component
    histfit(mv(real(residuals.tip(:,1,:)),imag(residuals.tip(:,1,:))),nbins)
    xlabel('Normalized Residual'); ylabel('Number of Data Points')
    title(['Tipper Tx: RMS = ',num2str(rms_comp(5))])

    subplot(3,3,9) %TY component
    histfit(mv(real(residuals.tip(:,2,:)),imag(residuals.tip(:,2,:))),nbins)
    xlabel('Normalized Residual'); ylabel('Number of Data Points')
    title(['Tipper Ty: RMS = ',num2str(rms_comp(6))])
    
else
    
    disp('Plotting impedance residuals')
    
    %Plot overall residual histograms
    subplot(3,3,[1 4 7])
    histfit([real(residuals.tip(:));imag(residuals.tip(:));real(residuals.imp(:));imag(residuals.imp(:))],nbins); hold on
    xlabel('Normalized Residual'); ylabel('Number of Data Points')
    title(['All Datapoints: RMS = ',num2str(s.rms)])

    subplot(3,3,2) %XX component
    histfit(mv(real(residuals.imp(:,1,:)),imag(residuals.imp(:,1,:))),nbins)
    xlabel('Normalized Residual'); ylabel('Number of Data Points')
    title(['Impedance XX: RMS = ',num2str(rms_comp(1))])

    subplot(3,3,3) %%XY component
    histfit(mv(real(residuals.imp(:,2,:)),imag(residuals.imp(:,2,:))),nbins)
    xlabel('Normalized Residual'); ylabel('Number of Data Points')
    title(['Impedance XY: RMS = ',num2str(rms_comp(2))])

    subplot(3,3,5) %YX component
    histfit(mv(real(residuals.imp(:,3,:)),imag(residuals.imp(:,3,:))),nbins)
    xlabel('Normalized Residual'); ylabel('Number of Data Points')
    title(['Impedance YX: RMS = ',num2str(rms_comp(3))])

    subplot(3,3,6) %%YY component
    histfit(mv(real(residuals.imp(:,4,:)),imag(residuals.imp(:,4,:))),nbins)
    xlabel('Normalized Residual'); ylabel('Number of Data Points')
    title(['Impedance YY: RMS = ',num2str(rms_comp(4))])

    subplot(3,3,8) %TX component
    histfit(mv(real(residuals.tip(:,1,:)),imag(residuals.tip(:,1,:))),nbins)
    xlabel('Normalized Residual'); ylabel('Number of Data Points')
    title(['Tipper Tx: RMS = ',num2str(rms_comp(5))])

    subplot(3,3,9) %TY component
    histfit(mv(real(residuals.tip(:,2,:)),imag(residuals.tip(:,2,:))),nbins)
    xlabel('Normalized Residual'); ylabel('Number of Data Points')
    title(['Tipper Ty: RMS = ',num2str(rms_comp(6))])

    
end



