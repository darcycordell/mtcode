function [new_err] = apply_errorfloor(tru_err,err_flr,denom)
% Function which takes a vector of real-valued errors and applies an error
% floor using some percentage of another quantity.
%
% Inputs:
%
%   tru_err: vector of errors
%   err_flr: some error floor (e.g. 0.05)
%   denom: the quantity which the errors are compared to (i.e. the
%   denominator)
%
%   The floor is applied such that if 
%
%           tru_err/denom < err_flr
%
%   then the new error is set as err_flr*denom
%
% In most contexts, this is used on complex impedance data, Z, with real-valued
% errors, Zerr. The denominator in most cases should be abs(|Zxy|*|Zyx|) although
% other variations are possible.

new_err = tru_err;
for flo=1:length(tru_err)
    if tru_err(flo)<err_flr*denom(flo)
        new_err(flo)=err_flr*denom(flo);
    end
end