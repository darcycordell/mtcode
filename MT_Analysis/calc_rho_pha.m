function [rho, pha, rhoerr, phaerr] = calc_rho_pha(Z,Zerr,T)
%Function which calculates apparent resistivity and phase from an impedance
%which is in SI units with e^{+iwt} convention.
%
% Usage: [rho, pha, rhoerr, phaerr] = calc_rho_pha(Z,Zerr,T)
%
%
% Inputs: Z and Zerr are nf x nr x ns impedances and impedance errors
%   T is a vector of nf x 1 periods
%
% Outputs: rho and pha are nf x nr x ns apparent resistivity and phases
%   rhoerr and phaerr are their respective nf x nr x ns errors.

mu = 4*pi*10^-7;
w = 2*pi./T;

sze = size(Z);

dZ = zeros(sze);
rhoerr = zeros(sze); rho = zeros(sze);
phaerr = zeros(sze); pha = zeros(sze);
for i = 1:sze(1)
    
    %The formula which converts impedance to apparent resistivity and phase
    rho(i,:,:) = (1/(w(i)*mu))*abs(Z(i,:,:)).^2;
    pha(i,:,:) = atan2(imag(Z(i,:,:)),real(Z(i,:,:)))*(180/pi);
    
    %The formulas which convert impedance error to error in apparent
    %resistivity and phase. Greg and Juliane wrote a PDF document with the
    %derivation.
    %Use real(Zerr)
    dZ(i,:,:) = real(real(Zerr(i,:,:))./sqrt(Z(i,:,:).*conj(Z(i,:,:))));
      
    rhoerr(i,:,:) = abs(2*rho(i,:,:).*(dZ(i,:,:)));
    %rhoerr(i,:,:) = rho(i,:,:).*((2*dZ(i,:,:))+(dZ(i,:,:)).^2);
    phaerr(i,:,:) = (180/pi)*dZ(i,:,:); %Phase error in degrees
    
end

