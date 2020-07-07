function convert_edi_spectra_imp
%
% Function which converts EDI files from Phoenix spectra format to impedance format
% Calls routine readEDI_spectra.m provided by Christina Walter
% All EDI files in folder are converted and placed in new folder edi_imp
% July 31 2016 : Compared to conversion in Winglink and results agree

% Still needed :  (a) error calculation
% (b) Check if results below are for remote reference or ocal processing, since
% both are possible.
%
% 

clear all; close all;

% X =1  Y = 2
cf = 1.25660E-03;  % 4*pi*1e-4

mkdir edi_imp

% List all *.edi files in folder
d=dir;
fff=length(d);
count=1;
for iio=1:fff
    [pathstr, name, ext] = fileparts(d(iio).name);
    if strcmp(ext,'.edi')
        edst{count}=d(iio).name;
        count=count+1;
    else
        disp([d(iio).name,' is not an edi file']);
    end
end

coord = []; pers = [];   nsta=length(edst);

% Loop over all EDI files
for is=1:nsta

  disp(['Reading ',edst{is}])
  [header,freq,imp_ten,tip_comp] = readEDI_spectra(edst{is});

  nfreq = length(freq);  
  lat = header.lat; long = header.lon; elev = header.elev;
  filename = header.name;
  
  % Reformat
  for ifreq = 1:nfreq
    Z(1,1,ifreq) =cf*complex(imp_ten(ifreq,1),imp_ten(ifreq,2));
    Z(1,2,ifreq) =cf*complex(imp_ten(ifreq,3),imp_ten(ifreq,4));
    Z(2,1,ifreq) =cf*complex(imp_ten(ifreq,5),imp_ten(ifreq,6));
    Z(2,2,ifreq) =cf*complex(imp_ten(ifreq,7),imp_ten(ifreq,8));
    Zvar(1,1,ifreq) = 0.01;   Zvar(1,2,ifreq) = 0.01;
    Zvar(2,1,ifreq) = 0.01;   Zvar(2,2,ifreq) = 0.01;
 
    T(1,ifreq) = complex(tip_comp(ifreq,1),tip_comp(ifreq,2));
    T(2,ifreq) = complex(tip_comp(ifreq,3),tip_comp(ifreq,4));
    Tvar(1,ifreq) = 0.01;     Tvar(2,ifreq) = 0.01;
 
    zrot(ifreq) = 0.;  trot(ifreq) = 0.; rhorot(ifreq) = 0.;
 
  end

  n22 = NaN(2,2,nfreq);

  % Write in impedance format
  write_edi_imp(Z,Zvar,n22,n22,n22,n22,T,Tvar,freq,lat,long,elev,zrot,rhorot,trot,['edi_imp/',edst{is}])

end