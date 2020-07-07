function [d]=rotate_d(d,rot)
%
% This function takes a standard data structure and rotates ALL data
% (impedance and tipper) by a single rotation angle. Positive rotation
% angle is clockwise and negative rotation angle is anti-clockwise.
%
% Usage: [d] = rotate_d(d,rot)
%
% Inputs: 
%
% "d" is standard MT data structure
% "rot" is the angle to rotate clockwise in degrees
%
% Note: The "d" data structure has impedances and tipper data which have
% already been rotated by d.zrot and d.trot. So for example, if d contains
% some impedance d.Z and d.zrot is 19 degrees, then d.Z has been rotated to
% 19 degrees. If you then do [d_new] = rotate_d(d,10), the data will be
% rotated by 10 degrees clockwise meaning that d_new will have impedances
% which have been rotated to 29 degrees. This code also updates the d_new.zrot
% and d_new.trot variables as well.

for is = 1:d.ns

    for ifreq = 1:d.nf
        
        d.Z(ifreq,:,is) = rotate_Z(d.Z(ifreq,:,is),rot);
        d.tip(ifreq,:,is) = rotate_tip(d.tip(ifreq,:,is),rot);
        
    end

end

[d.rho, d.pha, d.rhoerr, d.phaerr] = calc_rho_pha(d.Z,d.Zerr,d.T);

d.zrot = d.zrot + rot; % d.zrot should be nf x ns
d.trot = d.trot + rot;

% Leave errors unchanged. Fix later.
%dz = (U.^2*dz.^2).^.5;
%dz = abs(U*dz);
end

