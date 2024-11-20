function plot_fwd_1d(fwd)

T = 1./(fwd.freq_array);

subplot(2,2,1);
loglog(T,fwd.rho,'-k'); hold on
grid on
xlabel('Period (s)')
ylabel('App Rho (Ωm)')
set(gca,'FontSize',14)


subplot(2,2,3);
semilogx(T,fwd.phi,'-k'); hold on
grid on
xlabel('Period (s)')
ylabel('Phase (deg)')
set(gca,'FontSize',14)

subplot(2,2,2);
loglog(T,real(fwd.Z),'-k'); hold on
grid on
xlabel('Period (s)')
ylabel('Real Impedance (Ω)')
set(gca,'FontSize',14)

subplot(2,2,4);
loglog(T,imag(fwd.Z),'-k'); hold on
grid on
xlabel('Period (s)')
ylabel('Imag Impedance (Ω)')
set(gca,'FontSize',14)