function [worse,better,up,down,N]=plot_misfit_cross_plot(d,d1,d2,flag)
%Function to plot a cross plot of residuals from two different inversions
% This function is useful when comparing two forward responses to examine
% how different they are (e.g. when perturbing a model and comparing the
% original inversion respone to the perturbed forward response).
%
% Usage: [worse,better,up,down,N]=plot_misfit_cross_plot(d,d1,d2,flag)
%
% Inputs:
%       d: The observed data structure
%       d1: The first forward response data structure
%       d2: The second forward response data structure
%       flag: To show figures or not (1 = show, 0 = do not show)
%
% Outputs:
%       The number of data points which are worse, better, biased up,
%       biased down.
%       N: Total number of data points.

close all
[s1] = detailed_statistics(d,d1);
[s2] = detailed_statistics(d,d2);

d1.Zerr = d.Zerr;
[s3] = detailed_statistics(d1,d2);

%Make a vector residuals for each dataset
is = 1:d.ns;
ip = 1:d.nf; ir = 1:4;
a = mv(real(s1.residuals.imp(ip,ir,is)),imag(s1.residuals.imp(ip,ir,is))); %Residuals (d - d1)
b = mv(real(s2.residuals.imp(ip,ir,is)),imag(s2.residuals.imp(ip,ir,is))); %Residuals (d - d2)
c = mv(real(s3.residuals.imp(ip,ir,is)),imag(s3.residuals.imp(ip,ir,is))); %Residuals (d1 - d2)

a(isnan(a))=[];
b(isnan(b))=[];
c(isnan(c))=[];

if size(a) ~= size(b)
    disp('Error: number of data points must be equal for both datasets in order to plot cross-plot.')
    return
end

%avg = mean(c); %Usually approximately 0
%stddev = std(c); %Approximately the same as the stats_c.rms
xint = 1;%std(a);
yint = 1;

if flag
    set_figure_size(1);
    %Plot b vs c
    plot(a,b,'.k'); hold on
    plot(-10^6:10^6,-10^6:10^6,'--k') %Plot y=x line. Points near this line are points which began with low misfit and ended with high misfit
    plot((-10^6:10^6)-1,-10^6:10^6,'--k') %Plot y = 2x line. Points above this line had improved misfit because y = x+a where a = x.
    plot((-10^6:10^6)+1,-10^6:10^6,'--k') 
    plot(-10^6:10^6,-(-10^6:10^6),'--k') 

    %Plot the mean(c) line (usually close to y = 0) and the x=0 line
    plot([-10^6 10^6],[0 0],'-k','LineWidth',2) 
    plot([0 0],[-10^6 10^6],'-k','LineWidth',2)

    %Plot the y-interval. Here I use +/- 1. If points are outside these bounds,
    %it implies that they are larger than the error bars associated with those
    %points
    plot([-10^6 10^6],[yint yint],'--k')
    plot([-10^6 10^6],[-yint -yint],'--k')

    %Plot the x-interval. Here I use +/- 1. If points are outside these bounds, 
    %it implies they are larger than the error bars.
    plot([xint xint],[-10^6 10^6],'--r')
    plot([-xint -xint],[-10^6 10^6],'--r')

    fac = 6;
    axis([-fac*std(a) fac*std(a) -fac*std(b) fac*std(b)]); axis equal
    ylabel('Observed - Perturbed Residual (\Omega )')
    xlabel('Observed - Predicted Residual (\Omega )')
    
    plot(a,b,'.k'); hold on

end


k1 = 1; up_worse = [];
k2 = 1; down_worse = [];
k3 = 1; up_better = [];
k4 = 1; down_better = [];
for i = 1:length(a);
    
    if c(i)>1
        if b(i)>abs(a(i)) && b(i)>1
            up_worse(k1) = i;
            k1 = k1+1;
        end
        
        if b(i)<abs(a(i)) && a(i)<-1
            up_better(k3) = i;
            k3 = k3+1;
        end
    end
    
    if c(i)<-1
        
        if a(i)<abs(b(i)) && b(i)<-1
            down_worse(k2) = i;
            k2 = k2+1;
        end
        
        if a(i)>abs(b(i)) && a(i)>1
            down_better(k4) = i;
            k4 = k4+1;
        end
        
    end
    
end

worse = length(down_worse)+length(up_worse);
better = length(down_better)+length(up_better);
up = length(up_worse)+length(up_better);
down = length(down_worse)+length(down_better);
N = length(a);

if flag
    plot(a(up_worse),b(up_worse),'.b'); hold on
    plot(a(up_better),b(up_better),'.r'); hold on
    plot(a(down_worse),b(down_worse),'.b'); hold on
    plot(a(down_better),b(down_better),'.r'); hold on

    h(1) = plot(nan,nan,'.k'); hold on
    h(2) = plot(nan,nan,'.b');
    h(3) = plot(nan,nan,'.r');
    % h(4) = plot(nan,nan,'om');
    % h(5) = plot(nan,nan,'.g','MarkerSize',12);

    legend(h,'Residuals within Error','Worse','Better')
    title(['N = ',num2str(N),'. Improved: ',num2str(100*better/N),'%. Worse: ',num2str(100*worse/N),'%. Biased Up = ',num2str(100*up/N),'%. Biased Down = ',num2str(100*down/N),'%'])

%     print_figure(['statistics_compare'],['ks_test_cross_plot']); %Save figure
end

end % end function


