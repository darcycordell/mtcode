function [de] = cut_beta(d,betacutoff)
%
% Function which removes data above a given threshold of phase tensor beta
% skew
%
% Usage:
%       [de] = cut_beta(d,betacutoff)
%
% Inputs:
%           d = standard MT data structure
%
% Outputs:
%           de = standard MT data structure with edited points
%
%

count = 1;
for is = 1:d.ns  % Loop over stations

    [p] = calc_phase_tensor(d.Z(:,:,is)); %Calculate phase tensor parameters

    for ip = 1:d.nf
        if abs(p.beta(ip))>betacutoff
            K(count) = is;
            I(count) = ip;
            count = count+1;
            
        end
    end
    

end


for i = 1:length(K)
    d.rho(I(i),:,K(i)) = NaN;
    d.pha(I(i),:,K(i)) = NaN;
    d.Z(I(i),:,K(i)) = NaN+1i*NaN;
    d.rhoerr(I(i),:,K(i)) = NaN;
    d.phaerr(I(i),:,K(i)) = NaN;
    d.Zerr(I(i),:,K(i)) = NaN+1i*NaN;
end

de = d;