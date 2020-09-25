function dn = add_noise(d,flag)
%Function which adds random Gaussian noise and sets error bars for a given
%data set and writes out a new file with the noise added.
%
% Usage: dn = add_noise(d,flag)
%
% Inputs:
%   "d" is a data structure
%   flag is either 1 (to save as matfile d structure), 2 (to save outputs as 2-D), 
%           3 (to save as ModEM) or 0 (to not save outputs)
%       Default is to not save (flag = 0)
%
% Outputs:
%   "dn" is the data with noise output
%
% This function adds independent Gaussian noise to the real
% and imaginary components of the complex impedance. The "percent" Gaussian
% noise is taken to be either a percent of sqrt(|Zxy||Zyx|) for all
% components, or a percent of |Z_ij| for each ij component. In general, it
% is more consistent to use sqrt(|Zxy||Zyx|), although some cases (for
% example 1-D MT data) requires |Z_ij|.
%
%
% 

close all;

if ~exist('flag','var')
    flag = 0;
end

dn = d;

%Add NOISE level. This is the amount of "scatter" you want to add to the
%data (Gaussian noise)
noise_level={'0.05','0.02'};
gn=char(inputdlg({'Impedance noise level','Tipper noise level'},'Noise',1,noise_level));
GN = str2double(gn(1,:)); %Gaussian noise to be added to impedance data.
GN_tip = str2double(gn(2,:));

%GN = 0.02;
if isnan(GN)
    return
end


%%
%Add noise to impedance----------------------------------------------------
  
%Noise must be independent between real and imaginary components
alpha = randn(d.nf,4,d.ns); %random Gaussian matrix for real components
beta = randn(d.nf,4,d.ns); %random Gaussian matrix for imaginary components



option = menu('','Percent of |Z_xy| (DEFAULT)','Percent of |Z_ij|','Percent of sqrt(Zxy*Zyx)');
%option = 2;
if option == 3
%There are different ways to define the "noise percent". This one uses the
%percentage of sqrt(Z_xy*Z_yx) which makes the noise percentage consistent
%with the way inversions often define error floor (i.e. the r.m.s. between
%the noisy and noise-free data will be 1 if you apply an error floor in the
%normal way)
snr = GN*sqrt(abs(d.Z(:,2,:)).*abs(d.Z(:,3,:)));
snr(:,2,:) = snr(:,1,:);
snr(:,3,:) = snr(:,1,:);
snr(:,4,:) = snr(:,1,:);

dn.Z = real(d.Z)+alpha.*snr + 1i*(imag(d.Z)+beta.*snr);

elseif option == 2
%This method defines "noise percent" as a fraction of an individual
%impedance value where real and imaginary components are independent. To
%me, this seems most correct in that the variables are truly Gaussian and
%independent. However, this method makes it impossible to recover an rms of
%1.0 misfit between noisy and noise-free data when using ModEM error floors
%defined a different way. 
%
%****This method will give diagonal components which are not very noisy which may not be realistic.***
%
%
%snr = GN*d.Z; %signal-to-noise ratio
snr = GN*abs(d.Z); %signal to noise ratio
dn.Z = real(d.Z)+alpha.*snr + 1i*(imag(d.Z)+beta.*snr); %noisy impedance

else
%This method defines "noise percent" as a fraction of the off-diagonal
%components linked by E-fields. So Zxx and Zxy are defined relative to Zxy
%and Zyx and Zyy are defined relative to Zyx. This is the most complicated
%method but arguably gives the best results: diagonals are noisy (like real
%data), and error floors can still be used to recover an r.m.s. of 1.
    snr = GN*abs(d.Z);
    
    dn.Z = real(d.Z)+alpha.*snr + 1i*(imag(d.Z)+beta.*snr); %noisy impedance
    
    dn.Z(:,1,:) = real(d.Z(:,1,:))+alpha(:,1,:).*snr(:,2,:)+1i*(imag(d.Z(:,1,:))+beta(:,1,:).*snr(:,2,:));
    dn.Z(:,4,:) = real(d.Z(:,4,:))+alpha(:,4,:).*snr(:,3,:)+1i*(imag(d.Z(:,4,:))+beta(:,4,:).*snr(:,3,:));

end

%Define the error bars-----------------------------------------------------
dn.Zerr = zeros(size(dn.Z));
for i = [2 3]
    
    if option == 3
        dn.Zerr(:,i,:) = GN*sqrt(abs(d.Z(:,2,:)).*abs(d.Z(:,3,:)));
    
    elseif option == 2
        dn.Zerr(:,i,:) = GN*abs(d.Z(:,i,:)); %Errors as % of |Z_ij|
        
    else
        dn.Zerr(:,i,:) = GN*abs(d.Z(:,i,:)); %Errors as % of Z_xy and Z_yx
        
    end
    
    %Do not use this way:
    %dn.Zerr(:,i,:) = GN*d.Z(:,i,:); %Alternative method to define errors
    
end

for i = [1 4]
    
    if option == 3
        dn.Zerr(:,i,:) = GN*sqrt(abs(d.Z(:,2,:)).*abs(d.Z(:,3,:)));
        
    elseif option == 2
        dn.Zerr(:,i,:) = GN*abs(d.Z(:,i,:)); %Errors as % of |Z_ij|
        
    else
        
        dn.Zerr(:,1,:) = GN*abs(d.Z(:,2,:)); %Errors as % of Z_xy for Z_xx
        dn.Zerr(:,4,:) = GN*abs(d.Z(:,3,:)); %Errors as % of Z_yx for Z_yy
        
    end
    
    %Do not use this way:
    %dn.Zerr(:,i,:) = GN*d.Z(:,i,:); %Alternative method to define errors
    
end 

%Alternative way to define errors based on error propagation of real and
%imaginary errors where real error is GN*abs(real(d.Z)) and imaginary error
%is GN*abs(imag(d.Z)). This is based on the idea that the error is one
%standard deviation in a Gaussian sense.
%dn.Zerr = GN*(sqrt(real(d.Z).^4+imag(d.Z).^4)./sqrt(real(d.Z).^2+imag(d.Z).^2));

dn.Zerr = dn.Zerr+1i*dn.Zerr;

%plot_misfit_residual_histograms(dn,d)

[dn.rho, dn.pha, dn.rhoerr, dn.phaerr] = calc_rho_pha(dn.Z,dn.Zerr,dn.T);
%%

%Add noise to tipper------------------------------------------------------
%Noise must be independent between real and imaginary components
alpha = randn(d.nf,2,d.ns); %random Gaussian matrix for real components
beta = randn(d.nf,2,d.ns); %random Gaussian matrix for imaginary components

%snr = GN_tip*abs(d.tip); %signal to noise ratio using relative noise
snr = GN_tip*ones(size(d.tip)); %signal to noise ratio using absolute noise
dn.tip = real(d.tip)+alpha.*snr + 1i*(imag(d.tip)+beta.*snr); %noisy tipper
dn.tiperr = snr+1i*snr; %Errors as % of |T|


[stats_all] = detailed_statistics(dn,d);

disp(['RMS between noisy and noise-free datasets is: ',num2str(stats_all.rms)])

% plot_rho_pha(dn,25,u);
% plot_rho_pha(d,25,u)

if flag == 1
    %Save Matfile d structure
    datafile = [strtok(d.name,'.'),'_',num2str(GN),'_noise.mat'];
    d = dn;
    d.name = [d.name,'_noise'];
    save(datafile,'d');
elseif flag == 2
    prompt={'File Name','Data Mode (1 = TM; 2 = TE)'};
    dlg_title='File Name and Data Type';
    def={'Add_Noise.dat','1'};
    num_lines=1;
    dinp = inputdlg(prompt,dlg_title,num_lines,def);
    
    write_data_2DNLCG(dinp{1},dn,str2double(dinp{2}));
elseif flag == 3
    %Save ModEM File
    datafile = [strtok(d.name,'.'),'_',num2str(GN),'_noise.dat'];

    write_data_modem(datafile,dn)

else %If flag == 0 or something else, then do not save
    
end





