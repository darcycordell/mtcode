function [d] = load_data_edi_spectra(edi_file)

% Input: EDI file in Phoenix spectra format

% Output: d structure

% warning : rotation assumed zero and no errors calculated

% modified from original readEDI_spectra.m provided by Christina Walter

fid = fopen(edi_file);
key = ' ';
count=0;
while ~strcmp(key,'>END')
    st = textscan(fid,'%s',1);  
    key = cell2mat(st{1});
    if( size(key,2) > 6 && strcmp(key(1:7),'DATAID=') )     %  read name
        keylen = size(key,2);
        header.name = sscanf(key(8:keylen),'%s');
    end 
    if( size(key,2) > 6 && strcmp(key(1:4),'LAT=') )     %  read latitude
        keylen = size(key,2);
        ncoma = find(key == ':');
        deg = sscanf(key(5:ncoma(1)-1),'%d');
        minute = sscanf(key(ncoma(1)+1:ncoma(2)-1),'%d');
        sec = sscanf(key(ncoma(2)+1:keylen),'%f');
        header.lat = sign(deg)*(abs(deg)+minute/60+sec/3600);
    end 
    if( size(key,2) > 8 && strcmp(key(1:5),'LONG=') )     %  read long
        keylen = size(key,2);
        ncoma = find(key == ':');
        deg = sscanf(key(6:ncoma(1)-1),'%d');
        minute = sscanf(key(ncoma(1)+1:ncoma(2)-1),'%d');
        sec = sscanf(key(ncoma(2)+1:keylen),'%f');
        header.lon = sign(deg)*(abs(deg)+minute/60+sec/3600);
    end 
    if( size(key,2) > 5 && strcmp(key(1:5),'ELEV=') )     %  read long
        keylen = size(key,2);
        header.elev = sscanf(key(6:keylen),'%d');
    end 
    if( size(key,2) > 6 && strcmp(key(1:6),'NCHAN=') )     %  read number of channels
        keylen = size(key,2);
        header.nchan = sscanf(key(7:keylen),'%d');
    end 
    if( size(key,2) > 6 && strcmp(key(1:6),'NFREQ=') )     %  read number of frequencies
        keylen = size(key,2);
        header.nfreq = sscanf(key(7:keylen),'%i');
    end 
    if( size(key,2) > 5 && strcmp(key(1:5),'FREQ=') )     %  read frequencies
        keylen = size(key,2);
        count=count+1;
        freq(count,1) = sscanf(key(6:keylen),'%f');
    
        
%         nspec=cell2mat(textscan(fid,'%*s %*s %*s %*s %u',1)); % skip to number of spectra
        nspec=cell2mat(textscan(fid,'%*s %*s %*s %u',1)); % skip to number of spectra, edited for MTU-5C - exported from EMpower?

        
        xspec(count,:) = cell2mat(textscan(fid,'%f',nspec)); %read spectral density matrix in one row
        
    end   

    if( size(key,2) > 7 && strcmp(key(1:8),'>SPECTRA') ) 
        st = fgetl(fid);
        st = split(st);
        count=count+1;

        for sti = 1:length(st)
            if length(st{sti})>3 && strcmp(st{sti}(1:4),'FREQ')
                freq(count,1) = sscanf(st{sti}(6:end),'%f');
            end

            if length(st{sti})>2 && strcmp(st{sti}(1:3),'ROT')
                zrot(count,1) = sscanf(st{sti}(9:end),'%f');
            end
        end

        %nspec = sscanf(st{end},'%f');

        nspec = sscanf(cell2mat(join(regexp(st{end},'[0-9]','match'),'')),'%f');

        xspec(count,:) = cell2mat(textscan(fid,'%f',nspec));

    end


end

fclose(fid);

% check if data was read successfully
 E=exist('freq','var');
 
 if E==0
     
     freq=0;
     header=0;
     tip_comp=0;
     imp_ten=0;
     
 else
     % calculate impedance and tipper from spectra

     DET = ((xspec(:,36) - 1i.*xspec(:,6))  .* (xspec(:,44) - 1i.*xspec(:,14))  - (xspec(:,43) - 1i.*xspec(:,7))  .* (xspec(:,37) - 1i.*xspec(:,13)));
     ZXX = ((xspec(:,39) - 1i.*xspec(:,27)) .* (xspec(:,44) - 1i.*xspec(:,14))  - (xspec(:,46) - 1i.*xspec(:,28)) .* (xspec(:,37) - 1i.*xspec(:,13)))./DET;
     ZXY = ((xspec(:,46) - 1i.*xspec(:,28)) .* (xspec(:,36) - 1i.*xspec(:,6))   - (xspec(:,39) - 1i.*xspec(:,27)) .* (xspec(:,43) - 1i.*xspec(:,7)))./DET;
     ZYX = ((xspec(:,40) - 1i.*xspec(:,34)) .* (xspec(:,44) - 1i.*xspec(:,14))  - (xspec(:,47) - 1i.*xspec(:,35)) .* (xspec(:,37) - 1i.*xspec(:,13)))./DET;
     ZYY = ((xspec(:,47) - 1i.*xspec(:,35)) .* (xspec(:,36) - 1i.*xspec(:,6))   - (xspec(:,40) - 1i.*xspec(:,34)) .* (xspec(:,43) - 1i.*xspec(:,7)))./DET;

     
     TX = ((xspec(:,15) + 1i.*xspec(:,3))  .* xspec(:,9) - (xspec(:,16) + 1i.*xspec(:,10)) .* (xspec(:,8) + 1i.*xspec(:,2)))./DET;
     TY = ((xspec(:,16) + 1i.*xspec(:,10)) .* xspec(:,1) - (xspec(:,15) + 1i.*xspec(:,3))  .* (xspec(:,8) - 1i.*xspec(:,2)))./DET;
     
     
     
     imp_ten=[real(ZXX) imag(ZXX) real(ZXY) imag(ZXY) real(ZYX) imag(ZYX) real(ZYY) imag(ZYY)];

     tip_comp=[real(TX) imag(TX) real(TY) imag(TY)];

     
% Calculate and plot apparent resistivity and phase
      mu=4*pi*10^-7;
      
      
      rhoa_XY=mu./(2.*pi.*freq) .* (real(ZXY).^2 + imag(ZXY).^2).*10^6;
      rhoa_YX=mu./(2.*pi.*freq) .* (real(ZYX).^2 + imag(ZYX).^2).*10^6;
% 
      phase_XX= atan2d(imag(ZXX),real(ZXX)); 
      phase_XY= atan2d(imag(ZXY),real(ZXY));  
      phase_YX= atan2d(imag(ZYX),real(ZYX));
      phase_YY= atan2d(imag(ZYY),real(ZYY)); 
     
%      
% Plot results to verify     
      figure(2);
      subplot(4,1,1);
%      
      loglog(1./freq,rhoa_XY,'*r'); hold on;
      loglog(1./freq,rhoa_YX,'*b');
      grid on;
%       
      title(header.name);
      
      ylabel('Apparent resistivity [{\Omega}m]');
      xlim([0.003 10000]);
% 
%      
      subplot(4,1,2);
%      
      semilogx(1./freq,phase_XY,'*r');hold on;
      semilogx(1./freq,phase_YX,'*b');
      xlabel('Period [s]');
%      
      grid on;
%      
      %xlabel('Period [s]');
      ylabel('Phase [degrees]');
%      
      xlim([0.003 10000]);
      ylim([0 90]);
     
      ac=gca;
%     
      ac.YTick=([-90 -45 0 45 90]);
%     
      

      subplot(4,1,3);
      
      semilogx(1./freq,real(TX),'b*');hold on;
      semilogx(1./freq,-imag(TX),'r*');
  
      xlim([0.003 10000]);
      ylim([-1 1]);
      ylabel('Tx')
      
     
      subplot(4,1,4)
      semilogx(1./freq,real(TY),'b*'); hold on
      semilogx(1./freq,-imag(TY),'r*');
      
      xlim([0.003 10000]);
      ylim([-1 1]);
      ylabel('Ty')
      
      grid on
     
 end
 
 d.name = strrep(edi_file,'.edi','');
 d.niter = ' ';
 d.site = {d.name};
 d.T = 1./freq;
 d.f = freq;
 d.ns = 1;
 d.nf = numel(freq);
 d.nr = 4;
 d.responses = {'ZXX', 'ZXY', 'ZYX', 'ZYY', 'TX', 'TY'};
 d.Z = [ complex(imp_ten(:,1),imp_ten(:,2)) complex(imp_ten(:,3),imp_ten(:,4)) ...
 complex(imp_ten(:,5),imp_ten(:,6)) complex(imp_ten(:,7),imp_ten(:,8)) ] ./796; % put into SI units, standard for d structure
 d.Zerr = nan(size(d.Z)) + 1i.*nan(size(d.Z));
 
 [d.rho, d.pha, d.rhoerr, d.phaerr] = calc_rho_pha(d.Z,d.Zerr,d.T);

 d.tip = [complex(tip_comp(:,1),tip_comp(:,2)) complex(tip_comp(:,3),tip_comp(:,4))];
 d.tiperr = nan(size(d.tip)) + 1i.*nan(size(d.tip));
 d.zrot = zeros(d.nf, d.ns);
 d.trot = zeros(d.nf, d.ns);
 d.loc = [header.lat header.lon header.elev];
 d.x = 0;
 d.y = 0;
 d.z = header.elev;
 d.origin = d.loc;
 
 d = set_map_projection(d);