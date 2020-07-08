function [s] = detailed_statistics_2D(dobs,dpred,alpha)
%
% Function which calculates misfit statistics for almost every permutation
% possible given an observed data structure (nf x nr x ns) and predicted
% data structure (nf x nr x ns).
%
% Usage: [stats_all,rms_site,rms_freq,rms_comp,rms_sf,rms_fr,residuals] = detailed_statistics(dobs,dpred,alpha)
%
% Inputs:
% dobs - Standard data structure of observed data
% dpred - Standard data structure of predicted (modelled) data
% alpha = (OPTIONAL) the confidence level the user wishes to test 0<alpha<1
%
%
% Outputs:
% residuals - a nf x nr x ns matrix of normalized residuals for all data
%   points
%
% stats_all - a structure variable of statistics for the entire dataset
% (see statistics.m for details of the statistics structure).
%
% rms_site - a vector of rms values for each site
%
% rms_freq - a vector of rms values for each frequency
%
% rms_comp - a vector of rms values for each component (xx, xy, yx, yy, tx, ty)
%
% rms_sf - a ns x nf matrix of rms misfit values (sums over components)
%
% rms_fr - a nf x nr matrix of rms misfit values (sums over stations)
%
%


obs_vec = mv(dobs.rho,dobs.pha,real(dobs.tip),imag(dobs.tip));
pred_vec = mv(dpred.rho,dpred.pha,real(dpred.tip),imag(dpred.tip));
err_vec = mv(dobs.rhoerr,dobs.phaerr,real(dobs.tiperr),imag(dobs.tiperr));

s.residuals.rho=struct([]);
s.residuals.pha=struct([]);
s.residuals.tip=struct([]);

% residuals.imp = ((real(dobs.Z)-real(dpred.Z))./real(dobs.Zerr))+1i*((imag(dobs.Z)-imag(dpred.Z))./imag(dobs.Zerr));
s.residuals.rho = (dobs.rho - dpred.rho)./dobs.rhoerr;
s.residuals.pha = (dobs.pha - dpred.pha)./dobs.phaerr;
s.residuals.tip = ((real(dobs.tip)-real(dpred.tip))./real(dobs.tiperr))+1i*((imag(dobs.tip)-imag(dpred.tip))./imag(dobs.tiperr));

if exist('alpha','var')
    [stats] = statistics(obs_vec, pred_vec, err_vec,alpha);
    s.alpha = stats.alpha;
    s.chi_critical = stats.chi_critical;
    s.rms_critical = stats.rms_critical;
else
    [stats] = statistics(obs_vec, pred_vec, err_vec);
end
% put output variables into s structure
s.N = stats.N;
s.r = stats.r;
s.L2norm = stats.L2norm;
s.rms = stats.rms;
s.chisquare = stats.chisquare;

s.rms_site = zeros(dobs.ns,1); s.rms_sf = zeros(dobs.ns,dobs.nf);
for i = 1:dobs.ns

    obs_vec = mv(dobs.rho(:,:,i),dobs.pha(:,:,i),real(dobs.tip(:,:,i)),imag(dobs.tip(:,:,i)));
    pred_vec = mv(dpred.rho(:,:,i),dpred.pha(:,:,i),real(dpred.tip(:,:,i)),imag(dpred.tip(:,:,i)));
    err_vec = mv(dobs.rhoerr(:,:,i),dobs.phaerr(:,:,i),real(dobs.tiperr(:,:,i)),imag(dobs.tiperr(:,:,i)));

    [s.rms_site(i),s.N_site(i)] = calc_rms(obs_vec,pred_vec,err_vec);


    for j = 1:dobs.nf

        obs_vec = mv(dobs.rho(j,:,i),dobs.pha(j,:,i),real(dobs.tip(j,:,i)),imag(dobs.tip(j,:,i)));
        pred_vec = mv(dpred.rho(j,:,i),dpred.pha(j,:,i),real(dpred.tip(j,:,i)),imag(dpred.tip(j,:,i)));
        err_vec = mv(dobs.rhoerr(j,:,i),dobs.phaerr(j,:,i),real(dobs.tiperr(j,:,i)),imag(dobs.tiperr(j,:,i)));

        [rms] = calc_rms(obs_vec,pred_vec,err_vec);

        if isempty(rms)
            rms = NaN;
        end

        s.rms_sf(i,j) = rms;
    end

end

s.rms_freq = zeros(dobs.nf,1);
for i = 1:dobs.nf

    obs_vec = mv(dobs.rho(i,:,:),dobs.pha(i,:,:),real(dobs.tip(i,:,:)),imag(dobs.tip(i,:,:)));
    pred_vec = mv(dpred.rho(i,:,:),dpred.pha(i,:,:),real(dpred.tip(i,:,:)),imag(dpred.tip(i,:,:)));
    err_vec = mv(dobs.rhoerr(i,:,:),dobs.phaerr(i,:,:),real(dobs.tiperr(i,:,:)),imag(dobs.tiperr(i,:,:)));
    %aaa(i) = sum(~isnan(obs_vec));
    [s.rms_freq(i)] = calc_rms(obs_vec,pred_vec,err_vec);
end

%sum(aaa)

s.rms_comp = zeros(6,1); s.rms_fr = zeros(6,dobs.nf);
for i = 1:4 % XY and YX only non NaNs in 2D


    obs_vec = mv(dobs.rho(:,i,:),dobs.pha(:,i,:));
    pred_vec = mv(dpred.rho(:,i,:),dpred.pha(:,i,:));
    err_vec = mv(dobs.rhoerr(:,i,:),dobs.phaerr(:,i,:));

    [s.rms_comp(i)] = calc_rms(obs_vec,pred_vec,err_vec);

    for j = 1:dobs.nf

        obs_vec = mv(dobs.rho(j,i,:),dobs.pha(j,i,:));
        pred_vec = mv(dpred.rho(j,i,:),dpred.pha(j,i,:));
        err_vec = mv(dobs.rhoerr(j,i,:),dobs.phaerr(j,i,:));

        [rms] = calc_rms(obs_vec,pred_vec,err_vec);

        if isempty(rms)
            rms = NaN;
        end

        s.rms_fr(i,j) = rms;
    end


end

for i = 1:2
    [s.rms_comp(i+4)] = calc_rms([real(dobs.tip(:,i,:));imag(dobs.tip(:,i,:))],[real(dpred.tip(:,i,:));imag(dpred.tip(:,i,:))],[real(dobs.tiperr(:,i,:));imag(dobs.tiperr(:,i,:))]);

    for j = 1:dobs.nf

        obs_vec = mv(real(dobs.tip(j,i,:)),imag(dobs.tip(j,i,:)));
        pred_vec = mv(real(dpred.tip(j,i,:)),imag(dpred.tip(j,i,:)));
        err_vec = mv(real(dobs.tiperr(j,i,:)),imag(dobs.tiperr(j,i,:)));

        [rms] = calc_rms(obs_vec,pred_vec,err_vec);

        if isempty(rms)
            rms = NaN;
        end

        s.rms_fr(i+4,j) = rms;
    end

end

end