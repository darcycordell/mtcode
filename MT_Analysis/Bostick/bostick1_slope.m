function [rho_bostick_av,d_av, rho_berd_av,pha_berd_av] = bostick1(Z1,T)

% Basic Bostick transform to covert apparent resistivity and
% phase, to resistivity as a function of depth
% The method using here is calculating the slope of the apparent
% resistivity, so the phase is not used in this case at all

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
    %pha_berd_av_rad(ifreq) =  pha_berd_av(ifreq)/rad;  % phase in radians 

    % calculate the Bostick rho and corresponding depth
    d_av(ifreq)=sqrt(rho_berd_av(ifreq)/(w(ifreq)*mu));    %    in meter
    
end
ind=find(~isnan(rho_berd_av));

ninterp=101;
xxT= linspace(min(log10(T)),max(log10(T)),ninterp);
rho_temp=spline(log10(T(ind)),log10(rho_berd_av(ind)),xxT);
m_temp=gradient(rho_temp,xxT);

m=nan(1,nfreq);

for i=1:nfreq
    if ~isnan(rho_berd_av(i))
        difference=abs(xxT-log10(T(i)));
        ind2=find(difference==min(difference));
        if abs(m_temp(ind2))<1  %only when abs(m)<1, the bostick result is positive so we give a number to m. Or we keep nan for m
            m(i)=m_temp(ind2);
        end
    end
end  
rho_bostick_av=rho_berd_av.*(1+m)./(1-m);
d_av(isnan(rho_bostick_av))=NaN;


end