function [stats] = statistics(dobs,dpred,err,alpha)
%General function which calculates the chi-squared value for your dataset as well
%as the critical chi-squared value for your dataset given a confidence
%interval (alpha). Similarly for rms and L2-norm.
%
% Usage: [stats] = statistics(dobs,dpred,err,alpha)
%
%Inputs:
%
%dobs = vector or matrix of real-valued observed data
%dpred = vector or matrix of real-valued predicted data
%err = vector or matrix of real-valued errors (standard deviations)
%
%OPTIONAL:
%   alpha = the confidence level the user wishes to test 0<alpha<1
%
%Outputs:
%
%stats.chisquare = the chi-squared value for the given dataset
%stats.rms = the root-mean-square value for the given dataset
%stats.L2norm = the L2-norm for the dataset.
%
%OPTIONAL OUTPUTS when alpha is supplied:
%   stats.rms_critical = the critical RMS value.
%   stats.chi_critical = the critical chi-squared value. Chi-square values greater
%   than the critical value suggest that dpred is distinct from the noise of 
%   dobs with a confidence of alpha. In other words, there is only a (1-alpha)
%   chance that the two are the same given a random distribution.

dobs = dobs(:);
dpred = dpred(:);
err = err(:);

id = find(isnan(dobs)==1); %Find NaN indices. Note: in MT data, usually zeros 
        %represent null values or periods with no data. MT data impedances
        %cannot equal zero.
        
%Remove NaNs from vectors
dobs(id)=[];
dpred(id)=[];
err(id) = [];

%Remove zeros from vectors (in general, you should not remove zeros but,
%for MT data, it is unlikely you will have a data point that is exactly
%zero. Likely, if it is exactly zero, it is a non-data point).
id = find(dobs == 0);
if ~isempty(id)
    %disp(['Warning from ',mfilename,':'])
    %disp(['Observed data has ',num2str(numel(id)),' value(s) equal to zero. Check that this is correct! This might happen occasionally for tipper data.'])
end

% dobs(id)=[];
% dpred(id)=[];
% err(id)=[];

stats.N = length(dobs);

stats.r = (dobs-dpred)./err;

stats.L2norm = nansum(((dobs-dpred)./err).^2);

stats.rms = sqrt((1/stats.N)*nansum(((dobs-dpred)./err).^2));

stats.chisquare = nansum(((dobs-dpred)./err).^2);
%An optional way of computing chi-square (or rms) can be done using linear
%algebra
%       chisquare = (dobs - dpred).'*diag(1./err.^2)*(dobs - dpred)
%
%However, this method is about 50 times slower because it requires the
%creation of the NxN "data covariance matrix"

if exist('alpha','var') == 1
    
    stats.alpha = alpha;
    
    stats.chi_critical = chi2inv(alpha,stats.N);

    %RMS_from_chi = sqrt((1/N)*chisquare); %alternate defintion for RMS. Should
    %result in the same value.

    stats.rms_critical = sqrt((1/stats.N)*stats.chi_critical);

end


