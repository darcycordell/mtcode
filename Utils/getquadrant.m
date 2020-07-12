function [pha_edit]=getquadrant(pha)
% Function which puts phase data into the 0 - 90 quadrant
%
% Usage: [pha_edit] = getquadrant(pha)
%
% Inputs:
%       pha is a N x 1 vector of phase data
%
% Outputs:
%       pha_edit is an N x 1 vector of phase data with the phase data moved
%       to the first quadrant (0 - 90 degrees).
%
%

rad=pi/180;
if sin(nanmean(pha*rad))>0 && cos(nanmean(pha*rad))>0
    if nanmean(pha)<0
        pha_edit=pha+360;
    else%do nothing: phases are in correct quadrant!
        pha_edit=pha;
    end
elseif sin(nanmean(pha*rad))>0 && cos(nanmean(pha*rad))<0
    if nanmean(pha)<0
        pha_edit=pha+270;
    else
        pha_edit=pha-90;
    end
elseif sin(nanmean(pha*rad))<0 && cos(nanmean(pha*rad))<0
    if nanmean(pha)<0
        pha_edit=pha+180;
    else
        pha_edit=pha-180;
    end
elseif sin(nanmean(pha*rad))<0 && cos(nanmean(pha*rad))>0
    if nanmean(pha)<0
        pha_edit=pha+90;
    else
        pha_edit=pha-270;
    end
else
    if isnan(pha)
        pha_edit = NaN;
    else
        disp('Error: ***** There is a Serious Problem with you phases!!!*****')
    end
end

end