function convert_edi_units(filename)
%
% Function to convert EDI file from SI units to field units. 
%
% Usage: convert_edi_units(filename)
%
% Inputs: filename is the EDI file in SI units
%
% Outputs: Saves an EDI file in field units with filename_converted name
%
%
% All MTcode scripts assume an EDI file is in field units. If an EDI file is
% in SI units, then the EDI should be converted to field units before using
% "load_data_edi", "interpolate_edi", "MTplot" and other MTcode scripts.
%
% Impedance field units = E/B. E is in mV/km and B is in nT
%
% Impedance SI units = E/H. E is in V/m and H is in A/m
%
% To convert: Z_field = (Z_SI)/(mu*1000)
%
%



mu = 4*pi*10^-7;

[d] = load_data_edi(filename);

d.Z = d.Z/(mu*1000);
d.Zerr = d.Zerr/(mu*1000);
[d.rho,d.pha,d.rhoerr,d.phaerr]=calc_rho_pha(d.Z,d.Zerr,d.T); % calculate from the converted impedance

d.site = {[d.site{1},'_converted']};

disp(['Converting Site ',d.site{1},' from SI units to Field Units'])

if nanmedian(d.rho(:,2:3,:))<1 %If, after doing the conversion, the resistivies are still very small then something very strange is going on
    disp(['AFTER UNIT CONVERSION: Site ', filename(1:end-4),' STILL HAS VERY LOW RESISTIVITES (MEDIAN <1 OHM M). Please assess EDI file.'])
elseif nanmedian(d.rho(:,2:3,:))>10^6 %If after doing the conversion, the resistivites are huge then something else is very wrong. (Possibly the EDI file was already in field units?)
    disp(['AFTER UNIT CONVERSION: Site ', filename(1:end-4),' HAS VERY HIGH RESISTIVITES (MEDIAN >10^6 OHM M). Please assess EDI file.'])
else
    disp(['AFTER UNIT CONVERSION: Site ', filename(1:end-4),' has realistic resistivity values (Median = ',nanmedian(nanmedian(d.rho(:,2:3,:))),' Ohm m)'])
end

write_data_edi(d,1);




