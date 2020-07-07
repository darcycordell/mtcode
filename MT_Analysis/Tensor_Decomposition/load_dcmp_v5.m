%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load_dcmp - Program to read *.dcmp files from Strike program
%          - Modified from existing code in mj_pcolor_v4
%          - Output is a *.mat file with distortion parameters
% May 2008 - Ted Bertrand

% Version 5,   Feb 2011


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;


choice1 = menu('Choose periods to use','Load from file','Default values')
if choice1 == 1
    display('using user provided frequencies')
    [freq_file,freq_path] = uigetfile({'*.txt'},'pick  file containing frequencies');
    FDef = load( [freq_path freq_file], '-ascii');
else 
    display('using equally sampled periods')
    FDef=logspace(-4, 5, 72);  %IMPORTANT equally sampled frequencies
end

TDef = 1./FDef;

%Read in names of all files in working directory
files=ls; 
[numfiles maxfilename]=size(files);  
files=files(3:numfiles,:);             %eliminate hidden . and .. files
[numfiles maxfilename]=size(files);
numdcmp=0;                             %number of *.dcmp files

%Get the names of *.dcmp files (Note: Assume all *.dc* files are *.dcmp files)
for j=1:numfiles
    extension=find(files(j,:)=='.')+1;
    if extension>1
        files(j,extension:extension+1);
        if files(j,extension:extension+1)=='dc'  %keeps all files whose extension starts with 'dc'
            numdcmp=numdcmp+1;
            finst(numdcmp,:)=files(j,:);
        else
            disp([files(j,:),' is not a dcmp file']);
        end
    end
end
sname=finst;  %filenames of *.dcmp files

%Loop through each dcmp file and read in distortion parameters at each freq
[nsta,length]=size(sname);
for isite=1:nsta
    if exist(sname(isite,:),'file')
        fiddcmp=fopen(sname(isite,:),'r');
        disp(['File ',sname(isite,:),' opened']);
    else
        error(['File ',sname(isite,:),' not found!']);
        return
    end

    %Specify which parameters to read from DCMP files
    sparam = {'regional','shear','channelling','twist','rms','skew'};
    %===========================================
    % Extract parameters
    %===========================================
    for iextract=1:6
        line=1;
        while line ~= -1;
            line=fgetl(fiddcmp);
            if findstr(line,char(sparam{iextract}))
                nfreq(isite)=sscanf(line(1:4),'%f');
                for ifreq=1:nfreq(isite)
                    line=fgetl(fiddcmp);
                    temp2=sscanf(line,'%f');
                    per(ifreq)=temp2(1);
                    period(isite,ifreq)=temp2(1);
                    param(iextract,isite,ifreq)=temp2(2);
                end
                line=-1;
            end
        end
    end
    fclose(fiddcmp);
end
%period(find(period==0))=1e-40;

strike = squeeze(param(1,:,:));
shear = squeeze(param(2,:,:));
channelling = squeeze(param(3,:,:));
twist = squeeze(param(4,:,:));
misfit = squeeze(param(5,:,:));
skew = squeeze(param(6,:,:));

for j=1:numdcmp
strike_int(j,:)=interp1(log10(period(j,:)),(strike(j,:)),log10(TDef),'nearest');
shear_int(j,:)=interp1(log10(period(j,:)),(shear(j,:)),log10(TDef),'nearest');
channelling_int(j,:)=interp1(log10(period(j,:)),(channelling(j,:)),log10(TDef),'nearest');
twist_int(j,:)=interp1(log10(period(j,:)),(twist(j,:)),log10(TDef),'nearest');
misfit_int(j,:)=interp1(log10(period(j,:)),(misfit(j,:)),log10(TDef),'nearest');
skew_int(j,:)=interp1(log10(period(j,:)),(skew(j,:)),log10(TDef),'nearest');
end

%Uncomment to make plots and check interpolation
% for j=1:numdcmp
% semilogx(period(j,:),twist(j,:),'o');hold on;
% semilogx(TDef,twist_int(j,:),'*');hold on;
% end
% legend('original','interpolated')

%Write interpolated distortion parameter matrices to *.Mat file.
save([sname(1,1:3) '_dcmp'],'strike_int','shear_int','channelling_int','twist_int','misfit_int','skew_int');

