function [strike_max_split] = calc_pha_splits(d,is)
%%
% Function which calculates the strike angle which maximizes the splits
% between the TE and TM mode phase. This is done by looping through
% different rotation angles and calculating the phase split and then
% finding the maximum. This function is relatively slow. For a single
% station, it takes about 5 - 10 seconds.
%
% Usage: [strike_max_split] = calc_pha_splits(d,is)
%
% Inputs: 
%       "d" is standard MT data structure
%       "is" is station index to compute
%
% Outputs:
%       "strike_max_split" is the angle which the data must be rotated to
%       which achieves the greatest phase split

disp('Calculating Maximum Phase Split Angle')
nrot = 90; %Number of rotations to test

%Initialize matrices
dpha = nan(nrot,d.nf); rot = nan(nrot,1);
for irot = 1:nrot+1 %Loop through all rotation angles to test
  rot(irot) = (irot-1)*(180/nrot);
  [dr] = rotate_d(d,rot(irot));
  dpha(irot,:) = dr.pha(:,2,is)-(dr.pha(:,3,is)+180); %Phase split
end

%Determine the angle which contains the maximum phase split
ang_max = nan(d.nf,1);
ang_min = nan(d.nf,1);
for ifreq = 1:d.nf
    pdiff_max = -180;
    pdiff_min = 180;
    for irot = 1:nrot+1  %Loop through all rotations
        if dpha(irot,ifreq) >= pdiff_max %If the phase split for the irot rotation is larger than pdiff_max, then replace
            pdiff_max = dpha(irot,ifreq);
            ang_max(ifreq) = rot(irot);
        end

        if dpha(irot,ifreq) <= pdiff_min %If the phase split for irot rotation is smaller the pdfif_min, then replace
            pdiff_min = dpha(irot,ifreq);
            ang_min(ifreq) = rot(irot);
        end
    end
end

strike_max_split = ang_max;