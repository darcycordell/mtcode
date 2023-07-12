function d = append_d(d1,d2)
%Function to append one d structure dataset to another. Note that this only
%works if the two data structures how the identical frequency sets,
%responses, etc. It effectively appends new sites to an existing dataset.
%
% Usage: d = append_d(d1,d2)
%
%   Inputs:
%       d1: first dataset data structure
%       d2: second dataset data structure
%
%   Outputs:
%       d: appended dataset with [d1 d2]
%
%
d = d1;

d.Z = cat(3,d1.Z,d2.Z);
d.Zerr = cat(3,d1.Zerr,d2.Zerr);
d.tip = cat(3,d1.tip,d2.tip);
d.tiperr = cat(3,d1.tiperr,d2.tiperr);

d.site = cat(1,d1.site,d2.site);
d.loc = cat(1,d1.loc,d2.loc);

d.zrot = cat(2,d1.zrot,d2.zrot);
d.trot = cat(2,d1.trot,d2.trot);

d.rho = cat(3,d1.rho,d2.rho);
d.rhoerr = cat(3,d1.rhoerr,d2.rhoerr);
d.pha = cat(3,d1.pha,d2.pha);
d.phaerr = cat(3,d1.phaerr,d2.phaerr);

d.ns = size(d.Z,3);

d.name = 'Data_AB_BC_230310';

