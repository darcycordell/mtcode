function [rms,N] = calc_rms(dobs,dpred,err)
%General function which calculates only the root mean square error (rmse) value for a dataset
%
%       rms = sqrt((1/N)*sum(((dobs-dpred)/err)^2))
%
%
% Usage: [rms] = calc_rms(dobs,dpred,err)
%
%Inputs:
%
%dobs = vector or matrix of real-valued observed data
%dpred = vector or matrix of real-valued predicted data
%err = vector or matrix of real-valued errors (standard deviations)
%
%
%Outputs:
%
%rms = the root-mean-square value for the given dataset
%
%
% The function converts matrices to vectors and also automatically removes
% NaNs from the calculation.

dobs = dobs(:);
dpred = dpred(:);
err = err(:);

id = find(isnan(dobs)==1); %Find zero indices. Note: in MT data, usually zeros 
        %represent null values or periods with no data. MT data impedances
        %cannot equal zero.
        
%Remove NaNs from vectors
dobs(id)=[];
dpred(id)=[];
err(id) = [];

%Remove zeros from vectors (in general, you should not remove zeros but,
%for MT data, it is unlikely you will have a data point that is exactly
%zero. Likely, if it is exactly zero, it is a non-data point).

% Not always true for tipper data - the real or imaginary part could be
% zero. Commented out for now...
id = find(dobs == 0);
if ~isempty(id)
    %disp(['Warning from ',mfilename,':'])
    %disp(['Observed data has ',num2str(numel(id)),' value(s) equal to zero. Check that this is correct! This might happen occasionally for tipper data.'])
end
% 
% dobs(id)=[];
% dpred(id)=[];
% err(id)=[];

N = length(dobs);

rms = sqrt((1/N)*nansum(((dobs-dpred)./err).^2));

if isempty(rms)
    rms = NaN;
end



