function [d] = make_nan_data
%Function which makes a "standard" data structure full of NaN values at one
%frequency
%
% Usage: [d] = make_nan_data

%Header
d.name = 'None';
d.niter = '';
d.site = {'None'};

%Frequencies
d.T = 0;
d.f = 0;

%nf x nr x ns
d.ns = 1;
d.nf = 1;
d.nr = 4;

%Impedances
d.responses = {'None'};
d.Z = nan(d.nf,4,d.ns)+1i*nan(d.nf,4,d.ns);
d.Zerr = nan(d.nf,4,d.ns)+1i*nan(d.nf,4,d.ns);

%Apparent resistivity and phase
d.rho = nan(d.nf,4,d.ns);
d.rhoerr = nan(d.nf,4,d.ns);
d.pha = nan(d.nf,4,d.ns);
d.phaerr = nan(d.nf,4,d.ns);

%Tipper
d.tip = nan(d.nf,2,d.ns)+1i*nan(d.nf,2,d.ns);
d.tiperr = nan(d.nf,2,d.ns)+1i*nan(d.nf,2,d.ns);

%Rotations
d.zrot = [];
d.trot = [];

%Location and map info
d.loc = [0,0,0];
d.buffer = 0.1;
d.border = 1;
d.lim = [-0.1 0.1 -0.1 0.1];

%Model coordinate info
d.x = 0;
d.y = 0;
d.z = 0;
d.origin = [0,0];






