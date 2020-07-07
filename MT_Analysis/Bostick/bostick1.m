function [rho_bostick_av,d_av, rho_berd_av,pha_berd_av] = bostick1(Z1,T)

% Basic Bostick transform to covert apparent resistivity and
% phase, to resistivity as a function of depth
% The method using here is the simplized approximate way which uses the
% phase in the equation

% Written by Enci Wang, 2014
% Editted by MJU April 2014

nfreq = length(T);
mu = 4*pi*1e-7;rad = 180./pi; 
% Loop over MT frequencies
for ifreq = 1:nfreq  
    w(ifreq) = 2.*pi/T(ifreq);
end  

% Calculate the rho_bostick using the Berdichevsky average %%%%%%%%%%% 

for ifreq=1:nfreq
    Z_berd_av(ifreq) = (Z1(1,2,ifreq)- Z1(2,1,ifreq))/2;
    rho_berd_av(ifreq)= abs(Z_berd_av(ifreq)*Z_berd_av(ifreq))/(w(ifreq)*mu);
    pha_berd_av(ifreq)= 90-rad*atan2(real(Z_berd_av(ifreq)),imag(Z_berd_av(ifreq))); 
    pha_berd_av_rad(ifreq) =  pha_berd_av(ifreq)/rad;  % phase in radians
    if (pha_berd_av_rad(ifreq)>0) && (pha_berd_av_rad(ifreq)<pi/2)    %Only when phase is between 0 to 90 degree, the bostick result is positive
        % calculate the Bostick rho and corresponding depth
        d_av(ifreq)=sqrt(rho_berd_av(ifreq)/(w(ifreq)*mu)); 
        rho_bostick_av(ifreq)=rho_berd_av(ifreq)*(pi/(2*pha_berd_av_rad(ifreq))-1);
    else
        d_av(ifreq)=NaN;
        rho_bostick_av(ifreq)=NaN;
    end
end


end