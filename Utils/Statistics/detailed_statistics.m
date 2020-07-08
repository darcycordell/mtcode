function [s] = detailed_statistics(dobs,dpred,alpha)
%
% Function which calculates misfit statistics for almost every permutation
% possible given an observed data structure (nf x nr x ns) and predicted
% data structure (nf x nr x ns).
%
% Usage: [s] = detailed_statistics(dobs,dpred,alpha)
%
% Inputs:
% dobs - Standard data structure of observed data
% dpred - Standard data structure of predicted (modelled) data
% alpha = (OPTIONAL) the confidence level the user wishes to test 0<alpha<1
%
% Outputs:
% s - statistics structure with the following fields:
%
%     residuals - structure with two fields:
%          imp - a nf x nr x ns matrix of normalized impedance residuals
%          tip - a nf x nr x ns matrix of normalized tipper residuals
%
%     N - total number of data points
%
%     r - a column vector of normalized residuals of all non-NaN data points (i.e. length not
%     always = nf*nr*ns)
%
%     L2norm = the L2-norm for the entire dataset.
%
%     rms = the root-mean-square value for the entire dataset
%
%     chisquare = the chi-squared value for the entire dataset
%
%     rms_site - a vector of rms values for each site
%
%     rms_sf - a ns x nf matrix of rms misfit values (sums over components)
%
%     N_site - a 1 x ns vector containing the number of non-NaN data for each site
%
%     rms_freq - a vector of rms values for each frequency
%
%     rms_comp - a vector of rms values for each component (xx, xy, yx, yy, tx, ty)
%
%     rms_fr - a nf x nr matrix of rms misfit values (sums over stations)
%
%     ONLY IF INPUT CONTAINS ALPHA:
%
%     alpha - the input confidence interval
%    
%     chi_critical - the critical chi-squared value. Chi-square values greater
%     than the critical value suggest that dpred is distinct from the noise of 
%     dobs with a confidence of alpha. In other words, there is only a (1-alpha)
%     chance that the two are the same given a random distribution.
%
%     rms_critical - the critical RMS value.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obs_vec = mv(real(dobs.Z),imag(dobs.Z),real(dobs.tip),imag(dobs.tip));
pred_vec = mv(real(dpred.Z),imag(dpred.Z),real(dpred.tip),imag(dpred.tip));
err_vec = mv(real(dobs.Zerr),imag(dobs.Zerr),real(dobs.tiperr),imag(dobs.tiperr));

s.residuals.imp = ((real(dobs.Z)-real(dpred.Z))./real(dobs.Zerr))+1i*((imag(dobs.Z)-imag(dpred.Z))./imag(dobs.Zerr));
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

    obs_vec = mv(real(dobs.Z(:,:,i)),imag(dobs.Z(:,:,i)),real(dobs.tip(:,:,i)),imag(dobs.tip(:,:,i)));
    pred_vec = mv(real(dpred.Z(:,:,i)),imag(dpred.Z(:,:,i)),real(dpred.tip(:,:,i)),imag(dpred.tip(:,:,i)));
    err_vec = mv(real(dobs.Zerr(:,:,i)),imag(dobs.Zerr(:,:,i)),real(dobs.tiperr(:,:,i)),imag(dobs.tiperr(:,:,i)));

    [s.rms_site(i),s.N_site(i)] = calc_rms(obs_vec,pred_vec,err_vec);


    for j = 1:dobs.nf

        obs_vec = mv(real(dobs.Z(j,:,i)),imag(dobs.Z(j,:,i)),real(dobs.tip(j,:,i)),imag(dobs.tip(j,:,i)));
        pred_vec = mv(real(dpred.Z(j,:,i)),imag(dpred.Z(j,:,i)),real(dpred.tip(j,:,i)),imag(dpred.tip(j,:,i)));
        err_vec = mv(real(dobs.Zerr(j,:,i)),imag(dobs.Zerr(j,:,i)),real(dobs.tiperr(j,:,i)),imag(dobs.tiperr(j,:,i)));

        [rms] = calc_rms(obs_vec,pred_vec,err_vec);

        if isempty(rms)
            rms = NaN;
        end

        s.rms_sf(i,j) = rms;
    end

end

s.rms_freq = zeros(dobs.nf,1);
for i = 1:dobs.nf

    obs_vec = mv(real(dobs.Z(i,:,:)),imag(dobs.Z(i,:,:)),real(dobs.tip(i,:,:)),imag(dobs.tip(i,:,:)));
    pred_vec = mv(real(dpred.Z(i,:,:)),imag(dpred.Z(i,:,:)),real(dpred.tip(i,:,:)),imag(dpred.tip(i,:,:)));
    err_vec = mv(real(dobs.Zerr(i,:,:)),imag(dobs.Zerr(i,:,:)),real(dobs.tiperr(i,:,:)),imag(dobs.tiperr(i,:,:)));
    %aaa(i) = sum(~isnan(obs_vec));
    [s.rms_freq(i)] = calc_rms(obs_vec,pred_vec,err_vec);
end

%sum(aaa)

s.rms_comp = zeros(6,1); s.rms_fr = zeros(6,dobs.nf);
for i = 1:4


    obs_vec = mv(real(dobs.Z(:,i,:)),imag(dobs.Z(:,i,:)));
    pred_vec = mv(real(dpred.Z(:,i,:)),imag(dpred.Z(:,i,:)));
    err_vec = mv(real(dobs.Zerr(:,i,:)),imag(dobs.Zerr(:,i,:)));

    [s.rms_comp(i)] = calc_rms(obs_vec,pred_vec,err_vec);

    for j = 1:dobs.nf

        obs_vec = mv(real(dobs.Z(j,i,:)),imag(dobs.Z(j,i,:)));
        pred_vec = mv(real(dpred.Z(j,i,:)),imag(dpred.Z(j,i,:)));
        err_vec = mv(real(dobs.Zerr(j,i,:)),imag(dobs.Zerr(j,i,:)));

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