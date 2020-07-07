function [phase_tensor] = calc_phase_tensor(Z)
%
% Function which calculates the phase tensor and stores all the relevant
% parameters of the phase tensor in a structure.
%
% Usage: [phase_tensor] = calc_phase_tensor(Z)
%
% Inputs: "Z" is a (nf x 4) matrix of impedances for a single station
%
% Outputs: "phase_tensor" is a structure which contains the phase tensor,
% the maximum and minimum phase tensor ellipse axes, alpha angle, beta skew
% angle, ellipse eccentricity, geoelectric strike (i.e. alpha - beta), and
% the ellipse coordinates in x and y (for plotting purposes).
%
% See Caldwell et al., (2004) and Booker (2013) for a review of the phase
% tensor

rad = 180/pi; nrot = 30;

sze = size(Z);

X = real(Z);
Y = imag(Z);

%Initialize all matrices
pt = zeros(sze);
beta = zeros(sze(1),1);
alpha = zeros(size(beta));
phi1 = zeros(size(beta));
phi2 = zeros(size(beta));
phi3 = zeros(size(beta));
phimin = zeros(size(beta));
phimax = zeros(size(beta));

% Loop over all periods
for i=1:sze(1)
    
    % Eq. 13 of Caldwell et al. (2004)
    %PT(:,:,i) = ([X(i,1) X(i,2); X(i,3) X(i,4)]^-1)*[Y(i,1) Y(i,2); Y(i,3) Y(i,4)];
    
    % Alternatively, Eq. 15 of Caldwell et al. (2004)
    
    %Denominator
    detX = X(i,1)*X(i,4) - X(i,2)*X(i,3);
    
    %Calculate phase tensor using Eq. 15
    pt(i,1) = (X(i,4)*Y(i,1) - X(i,2)*Y(i,3))/detX;
    pt(i,2) = (X(i,4)*Y(i,2) - X(i,2)*Y(i,4))/detX;
    pt(i,3) = (X(i,1)*Y(i,3) - X(i,3)*Y(i,1))/detX;
    pt(i,4) = (X(i,1)*Y(i,4) - X(i,3)*Y(i,2))/detX;
    
    %Singular Value Decomposition gives same answer as above
    %[u,s,v] = svd([pt(i,1) pt(i,2); pt(i,3) pt(i,4)])

    
    % Eq. 19 of Caldwell et al. (2004)
    beta(i) = rad*1/2*atan2((pt(i,2)-pt(i,3)),(pt(i,1)+pt(i,4)));
    
    % Eq. 22 of Caldwell et al. (2004)
    alpha(i) = rad*1/2*atan2((pt(i,2)+pt(i,3)),(pt(i,1)-pt(i,4)));
    
    % Eq. A5 to A7 of Caldwell et al. (2004)
    phi1(i) = trace([pt(i,1) pt(i,2); pt(i,3) pt(i,4)])/2;
    phi2(i) = sqrt(det([pt(i,1) pt(i,2); pt(i,3) pt(i,4)]));
    phi3(i) = (pt(i,2)-pt(i,3))/2;
    
    %Compute the x and y coordinates for the phase tensor ellipse
    for irot=1:nrot
        rot_ang = (irot-1)*360/(nrot-1);
        v = [cosd(rot_ang),sind(rot_ang)];
        plot_pt = [pt(i,1) pt(i,2); pt(i,3) pt(i,4)]*v';
        x(irot,i) = plot_pt(1);    y(irot,i) = plot_pt(2);  
    end 
    
    % Find maximum value of phase ellipse at this frequency
    m1 = max(abs(x(:,i))); m2 = max(abs(y(:,i)));
    max_pt(i) = 1.1*max([m1,m2]);
    
end

% Eq. A8 and A9 of Caldwell et al. (2004)
phimin = rad*atan(sqrt(phi1.^2+phi3.^2)-sqrt(phi1.^2+phi3.^2-phi2.^2));
phimax = rad*atan(sqrt(phi1.^2+phi3.^2)+sqrt(phi1.^2+phi3.^2-phi2.^2));

e = sqrt(1-(phimin.^2./phimax.^2)); %Eccentricity of ellipse
strike = alpha - beta; %geoelectric strike

phase_tensor.pt = pt;
phase_tensor.phimin = phimin;
phase_tensor.phimax = phimax;
phase_tensor.alpha = alpha;
phase_tensor.beta = beta;
phase_tensor.e = e;
phase_tensor.strike = strike;
phase_tensor.x = x;
phase_tensor.y = y;
phase_tensor.max_pt = max_pt;
