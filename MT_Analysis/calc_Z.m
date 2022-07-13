function [Z,Zerr] = calc_Z(rho,rhoerr,pha,phaerr,T)
%
%Function which calculates SI unit impedance from apparent resistivity and
%phases assuming an e^{+iwt} convention
%
% [Z,Zerr] = calc_Z(rho,rhoerr,pha,phaerr,T)
%
%
% Whether you use apparent resistivity error or phase error to calculate
% the error in Z gives the same results. Here, we use apparent resistivity
% error.

mu = 4*pi*10^-7;
w = 2*pi./T;

sze = size(rho);

Zerr = zeros(sze); Z = zeros(sze);
for i = 1:sze(1)
    
    %The formula which converts apparent resistivity to impedance
        %NOTE: This uses a e^(+iwt) convention.
    Z(i,:,:) = sqrt(rho(i,:,:).*w(i)*mu).*exp(1i*(pha(i,:,:)*(pi/180)));
        
    %The formulas which convert impedance error to error in apparent
    %resistivity and phase and vice versa. Greg and Juliane wrote a PDF document with the
    %derivation.
    Zerr(i,:,:) = 0.5*(rhoerr(i,:,:)./rho(i,:,:)).*abs(Z(i,:,:));
    
    %This is an alternative definition using the phase error (gives
    %identical results as above)
    Zerr(i,:,:) = (pi/180)*phaerr(i,:,:).*abs(Z(i,:,:));
    
end

Zerr = Zerr + 1i*Zerr;
