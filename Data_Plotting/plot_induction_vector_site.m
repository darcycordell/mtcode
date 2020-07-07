function plot_induction_vector_site(d,is)
%%
% Function which plots the real and imaginary induction vectors at a given
% site as a function of frequency
%
% Usage: plot_induction_vector_site(d,is)
%
% "d" is data input in standard format
% "is" is the station index to plot
%
% The plots are in subplots 233 and 236 because this code is paired with
% plot_tipper

%Find real and imaginary tipper components
tzxr = real(d.tip(:,2,is));
tzxi = -imag(d.tip(:,2,is));
tzyr= real(d.tip(:,1,is));
tzyi= -imag(d.tip(:,1,is));

%Real induction vector
subplot(2,3,3)
for ifreq=1:2:d.nf %skip every second frequency
    
    x_orig = log10(d.T(ifreq));   y_orig = 0;
    % Plot induction vector as a line
    xplot =[x_orig,x_orig+2*tzyr(ifreq)]; %x coordinate
    yplot =[y_orig,y_orig+2*tzxr(ifreq)]; %y coordinate
    plot (xplot, yplot,'b-');  hold on;      
    plot (xplot(2),yplot(2),'b.')

end  
title('Real Induction Vector')
axis([min(log10(d.T)) max(log10(d.T)) -2 2])
xlabel('log10(Period (s))')
ylabel('Real IV')
grid on

%Imaginary induction vector
subplot(2,3,6)
plot([min(log10(d.T)) max(log10(d.T))], [0 0],':k'); hold on
for ifreq=1:2:d.nf 
    x_orig = log10(d.T(ifreq));   y_orig = 0;
    % Plot as line
    xplot =[x_orig,x_orig+2*tzyi(ifreq)]; 
    yplot =[y_orig,y_orig+2*tzxi(ifreq)]; 
    plot (xplot, yplot,'r-');  hold on;      
    plot (xplot(2),yplot(2),'r.')
end  
title('Imaginary Induction Vector')
axis([min(log10(d.T)) max(log10(d.T)) -2 2])
xlabel('log10(Period (s))')
ylabel('Imag IV')
grid on
