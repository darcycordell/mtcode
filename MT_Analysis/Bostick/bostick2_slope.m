
function [good,rho_bostick_av,d_av,rho_berd_av,pha_berd_av,rho_bostick_max,rho_bostick_min,d_max,max_dir,min_dir] = bostick2(Z1,Z1var,T,d_interp)

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

%%%%%%%%%%%%%%%%%%% Calculate the rho_bostick using the average %%%%%%%%%%% 

[rho_bostick_av,d_av, rho_berd_av,pha_berd_av] = bostick1(Z1,T);

%%%%%%%%%% Calculate the rho_bostick of all angles (1-181 degree) %%%%%%%%%

if nargin==3
    %use the scale of the depth from Bostick_av to define the depth of interpolation
    logdepth=log10(d_av);
    d_interp=min(min(logdepth)):0.01:max(max(logdepth));
end


% figure(33)
for irot=1:181
    
    %rotate impedances
    [Z1_rot,Z1var_rot] = rot_z_dr(Z1,Z1var,nfreq,irot);
    
    for ifreq=1:nfreq
        rho_xy(ifreq,irot)= abs(Z1_rot(1,2,ifreq)*Z1_rot(1,2,ifreq))/(w(ifreq)*mu);
%         pha_xy(ifreq,irot)= 90-rad*atan2(real(Z1_rot(1,2,ifreq)),imag(Z1_rot(1,2,ifreq))); 
%         pha_xy_rad(ifreq,irot) =  pha_xy(ifreq,irot)/rad;  % phase in radians 
        
        % calculate the Bostick rho and corresponding depth
        d(ifreq,irot)=sqrt(rho_xy(ifreq,irot)/(w(ifreq)*mu)); 
%         rho_bostick(ifreq,irot)=rho_xy(ifreq,irot)*(pi/(2*pha_xy_rad(ifreq,irot))-1);
    end
    
    ind=find(~isnan(rho_xy(:,irot)));
    
    ninterp=101;
    xxT= linspace(min(log10(T)),max(log10(T)),ninterp);
    rho_temp=spline(log10(T(ind)),log10(rho_xy(ind,irot)),xxT);
    m_temp=gradient(rho_temp,xxT);
    
    m(:,irot)=nan(nfreq,1);
    for i=1:nfreq
        if ~isnan(rho_xy(i,irot))
            difference=abs(xxT-log10(T(i)));
            ind2=find(difference==min(difference),1,'first');
            if abs(m_temp(ind2))<1  %only when abs(m)<1, the bostick result is positive so we give a number to m. Or we keep nan for m
                m(i,irot)=m_temp(ind2);
            end
        end
    end    
        
    rho_bostick_rot(:,irot)=rho_xy(:,irot).*(1+m(:,irot))./(1-m(:,irot));
    d(isnan(rho_bostick_rot))=NaN;
    clear rho_temp m_temp difference ind ind2;
    
    % Interpolate to uniform depth using two temporary variable x,y

      x = log10(d(:,irot));       x(isnan(x))=[]; % Remove the NaN values from the vector
      y = log10(rho_bostick_rot(:,irot));    y(isnan(y))=[]; % Remove the NaN values from the vector
  if length(x)>5    
      log_rho_bostick_rot(:,irot)=interp1(x,y,d_interp,'cubic',NaN);
      good=1;
  else
      good=0;
      rho_bostick_max=nan(size(d_interp));
      rho_bostick_min=nan(size(d_interp));
      d_max=nan(size(d_interp));
      max_dir=nan(size(d_interp));
      min_dir=nan(size(d_interp));
      return;
      
  end
    
%     plot(log_rho_bostick_rot,d_max,'-');hold on 
    
end

% pause; close(33)

d_max=10.^(d_interp);
rho_bostick_rot=10.^log_rho_bostick_rot;

%%%%%%%%%%%%%%%%% Select the maximum resistivity for each depth %%%%%%%%%%%
save bostick;
for idepth=1:length(d_max)
    ind=find(~isnan(rho_bostick_rot(idepth,:)));
    direction=1:181;    newdir=direction(ind);
    new_rho_bostick_rot=rho_bostick_rot(idepth,ind);
    if ~isempty(ind)
%       max_ind=find(new_rho_bostick_rot==max(new_rho_bostick_rot));
      max_ind=find(new_rho_bostick_rot==max(new_rho_bostick_rot),1,'first');
      rho_bostick_max(idepth)=new_rho_bostick_rot(max_ind);
      max_dir(idepth)=newdir(max_ind);
%       min_ind=find(new_rho_bostick_rot==min(new_rho_bostick_rot));
      min_ind=find(new_rho_bostick_rot==min(new_rho_bostick_rot),1,'first');
      rho_bostick_min(idepth)=new_rho_bostick_rot(min_ind);
      min_dir(idepth)=newdir(min_ind);
    else
      max_dir(idepth) = NaN;
      rho_bostick_max(idepth) = NaN;
      min_dir(idepth) = NaN;
      rho_bostick_min(idepth)=NaN;
    end

end
end