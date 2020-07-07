close all; clear all;
curdir = pwd;
data_menu = menu('','Synthetic Data','EDI','Text File');

if data_menu == 1
    
    %USER INPUTS FOR SYNTHETIC DATA
    datatype = 'SYNTHETIC';
    prompt={'Number of Frequencies','Minimum Frequency','Maximum Frequency','Gaussian Error to Add','Depth to Top of Model Layers','Layer Resistivites'};
    dlg_title='Synthetic Data Parameters';
    def={num2str(80),num2str(0.001),num2str(1000),num2str(0.05),'[0 150 250 5000 5500]','[100 1 50 1 10]'};
    %def={num2str(50),num2str(0.001),num2str(1000),num2str(0.1),'[0 2000 10000]','[100 1 1000]'};
    num_lines=1;
    dinp = inputdlg(prompt,dlg_title,num_lines,def);

    if isempty(dinp)
        return
    end

    num_freq = str2double(dinp{1});
    min_freq = str2double(dinp{2});
    max_freq = str2double(dinp{3});
    err = str2double(dinp{4}); %Gaussian error to add to data
    model_depth = str2num(dinp{5}); %Depth to top of each model layer
    model_res = str2num(dinp{6}); %Resistivity of each layer

    % FROM MT1D SYNTHETIC
    f = 10.^linspace(log10(min_freq),log10(max_freq),num_freq).';
    nd = length(f);
    T = 1./f;
    lat = 0;
    lon = 0;
    elev = 0;
    rot = 0;

    [fwd]=calc_fwd_1d(model_depth,model_res,f',err);

    Z = fwd.Z.';
    rhoa = fwd.rho;
    pha = fwd.phi;

    dZ = err*abs(Z);
    
elseif data_menu == 2
    %LOAD EDI
    [edifile,edipath] = uigetfile({'*.edi'},'Pick EDI file');

    if edifile == 0
        return
    end

    cd(edipath)
    DS = load_data_edi(edifile);
    lat = DS.loc(1);
    lon = DS.loc(2);
    elev = DS.loc(3);
    rot = median(DS.zrot);
    f = DS.f; %Frequencies (Hz)
    T = 1./f; %Periods (s)
    nd = length(f);
    cd(curdir)

    mode_menu = menu('Type of 1-D Data to Use:','XY mode','YX mode','Determinant Average (default)');

    if mode_menu == 1 %Fit the TE mode
        datatype = 'XY';
        rhoa = DS.rho(:,2)'; %apparent resistivity data
        Z = DS.Z(:,2); %impedance
        dZ = real(DS.Zerr(:,2)); %impedance error
        pha = DS.pha(:,2)'; %phase data
    elseif mode_menu == 2 %Fit the TM mode
        datatype = 'YX';
        rhoa = DS.rho(:,3)'; %swap modes
        pha = DS.pha(:,3)'; %swap modes
        %Depending on processing and EDI formatting, it is sometimes necessary
        %to add +180 to the TM EDI phases
        if mean(pha)>0 && mean(pha)<90
            %do nothing
        else
            pha = pha+180;
        end
        Z = DS.Z(:,3)*exp(1i*-pi); %Rotates yx impedance by 180 degree to make the TM mode "look" like the 1-D TE mode
        dZ = real(DS.Zerr(:,3)); %Not sure if I should rotate TM errors...
    else %Fit the determinant impedance
        datatype = 'DETERMINANT';
        %Determinant impedance as defined by Rung-Arunwan et al., 2017
        Z = (sqrt((DS.Z(:,1).^2+DS.Z(:,2).^2+DS.Z(:,3).^2+DS.Z(:,4).^2)./2));

        %The error is propagated using standard error propagation rules
        % A^2 --> dA^2
        % sqrt(A) --> sqrt(dA)
        % a+b+c --> sqrt(da^2 + db^2 + dc^2)
        dZ = real(sqrt(0.5*sqrt(DS.Zerr(:,1).^4+DS.Zerr(:,2).^4+DS.Zerr(:,3).^4+DS.Zerr(:,4).^4)));

        rhoa = ((1./(2*pi*DS.f'*4*pi*10^-7)).*abs(Z.').^2); %Apparent resistivity formula
        pha = (atan2(imag(Z),real(Z))*(180/pi)).'; %Phase formula

    end

    dZ_orig = dZ;

elseif data_menu == 3
    %LOAD TEXT FILE
    [datfile,datpath] = uigetfile({'*'},'Pick text file');
        
    if datfile == 0
        return
    end

    cd(datpath)
    [~,~,ext] = fileparts(datfile);

    if strcmp(ext,'.dat') || strcmp(ext,'.txt') || strcmp(ext,'.csv') || strcmp(ext,'.xls')
        d = table2array(readtable(datfile));
    elseif strcmp(ext,'.resp')       
        fid = fopen(datfile);
        d = fscanf(fid,'%e %e %e %e %e',[5 Inf])';
        fclose(fid);
    else
        error('Error: File type not supported. Must be .dat, .txt, .csv, .xls, or .resp')
    end

    cd(curdir)

    lat = 0;
    lon = 0;
    elev = 0;
    rot = 0;

    datatype = 'TXT';
    f = d(:,1); %Frequencies (Hz)
    T = 1./f; %Periods (s)
    nd = length(f);
    rhoa = d(:,2)'; %apparent resistivity data
    pha = d(:,3)'; %phase data
    rhoerr = d(:,4); %apparent resistivity error
    phaerr = d(:,5); %phase error

    [Z,dZ] = calc_Z(rhoa',rhoerr,pha',phaerr,T);
    
    dZ = real(dZ);
    
end

prompt={'D+ Smoothing Error Floor'};
dlg_title='D+ Smoothing';
def={'5'};
num_lines=1;
dinp = inputdlg(prompt,dlg_title,num_lines,def);
dplus_smoothing = str2double(dinp{1});


%-------------------------------------------------------------------------
d = dplus(Z,dZ+1i*dZ,T,dplus_smoothing);

for ifreq=1:length(T)
    if abs(abs((dZ(ifreq)))/abs(Z(ifreq))) < dplus_smoothing/100
        dZ(ifreq)=(dplus_smoothing/100)*abs(Z(ifreq));
        dZ(ifreq) = dZ(ifreq) + 1i*dZ(ifreq);
    end
end
[~,~,rhoerr,phaerr] = calc_rho_pha(Z,dZ,T);

set_figure_size(1);
subplot(2,2,1);
logerrorbar(T,rhoa',rhoerr,'.k','-k'); hold on
loglog(T,d.rho,'-k');
axis([min(T) max(T) 1 1000])
grid on
xlabel('Period (s)')
ylabel('\rho _a (\Omega m)')

subplot(2,2,3);
errorbar(T,pha',phaerr,'.k'); hold on
loglog(T,d.pha,'-k');
set(gca,'XScale','log');
axis([min(T) max(T) 0 90]);
grid on
xlabel('Period (s)');
ylabel('Phase (deg)')

subplot(2,2,2)
stairs([d.lambda(1) d.lambda 10^9],[10^-9 d.a_step d.a_step(end)],'k-')
set(gca,'XScale','log'); grid on
set(gca,'YScale','log');
axis([min(d.lambda)/10 max(d.lambda)*10 min(d.a_step)/10 max(d.a_step)*10])
xlabel('\lambda');
ylabel('Spectral Function a(\lambda)')
title(['R.M.S. = ',num2str(d.rms)])

subplot(2,2,4)
for k=1:length(d.z)
    loglog([10^-16 d.tau(k)],[d.z(k)/1000 d.z(k)/1000],'-k','LineWidth',1);hold on
end
axis([10^0 10^4 10^-1 10^2]); axis ij
set(gca,'Xtick',[10^0 10^1 10^2 10^3 10^4])
Xtick=[10^0 10^1 10^2 10^3 10^4];
XTickLabels = cellstr(num2str(round(log10(Xtick(:))), '10^%d'));
xlabel('Conductance (Siemens)')
ylabel('Depth (km)');
grid on
