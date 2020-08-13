function T = frequency_analysis
%Function which takes a set of EDI files in the current directory and determines 
%the frequency content of those EDIs. 
%It outputs two figures: 
%   1) A pcolor plot of frequencies vs. sites where black 
%       signifies data is present and white signifies no data is present
%   2) A histogram of all frequencies present.
%
% Usage: frequency_analysis
%
% Inputs: None
%
% Outputs: Outputs the frequency list
%
% Option to look at the raw frequencies of the EDIs or interpolate and then
% look at the frequencies. Note that the histograms for both raw and
% interpolated data sets should look roughly the same!
%
%

clear all; close all; 
curdir = pwd;
%Option to determine raw frequencies or interpolate
iopt = menu('Choose option','Determine Raw Frequencies Present','Interpolate EDIs','Matfile (d structure)');

if iopt == 1 %Option to load all EDI raw data------------------------------

    %Find all the EDI filenames in the current directory
    directory=dir;
    fff=length(directory);
    count=1;
    for iio=1:fff
        [~, ~, ext] = fileparts(directory(iio).name);
        if strcmp(ext,'.edi')
            edst{count}=directory(iio).name;
            count=count+1;
        else
            disp([directory(iio).name,' is not an edi file']);
        end 
    end

    nsta=length(edst); %nsta is the number of edi files in the directory
    
    %Loop over EDI files
    T = [];
    for iyu=1:nsta
        % Read in EDI file
        disp(['Reading ',edst{iyu}])
        [d(iyu)] = load_data_edi(edst{iyu});
        
        %Add all the frequencies to a running vector
        T = [T; d(iyu).T];
        site{iyu} = d(iyu).site;

    end
    
    T_all = T; %All frequencies present
    
    tol = 10^-5; %Logarithmic tolerance to determine if a frequency is "unique"
    T = 10.^(uniquetol(log10(T),tol)); %Determine unique frequency within tolerance
    uniqueT = T;
    
    %Determine which sites have which frequencies
    freq = zeros(length(T),length(d));
    for is = 1:length(d)
        for ifreq = 1:length(T)
            %If a particular site has a particular frequency then the
            %difference is less than the tolerance
            [freqa] = nearestpoint(log10(T(ifreq)),log10(d(is).T));
            [inda] = min(abs(log10(T(ifreq))-log10(d(is).T)))<=tol;
            if inda
                if ~all(isnan(d(is).Z(freqa,:,1)))
                    freq(ifreq,is) = inda; %Add this to the frequency matrix
                end
            end
        end
    end
    
    ns = length(d);
    
elseif iopt == 2 %Option to Interpolate EDIs-------------------------------
    
    d = interpolate_edi; %Interpolate EDI
    
elseif iopt == 3
    
    [matfile,matpath] = uigetfile({'*.mat'},'Pick Mat File');
    
    cd(matpath);
    load(matfile);
    cd(curdir);
    
else
    return
    
end

if iopt == 2 || iopt == 3
        
    %Determine which sites have which frequencies
    ns = d.ns;
    freq = NaN(d.nf,d.ns); T=[];
    for is = 1:d.ns
        for ifreq = 1:d.nf
            %If no data is present, then vals are NaN
            %[inda] = ~(all(isnan(d.Z(ifreq,:,is))) &&all(isnan(d.tip(ifreq,:,is)))); Tipper included 
            [inda] = ~(all(isnan(d.Z(ifreq,:,is))));
            freq(ifreq,is) = inda;
            if inda
                T = [T; d.T(ifreq)];
            end
        end
    end
    
    T_all = T;
    site = d.site;
    uniqueT = unique(T);
    T = d.T;
    
end


freq(:,ns+1) = 1;

[X,Y] = meshgrid(1:ns+1,log10(1./T));
  
%Plot a pcolor figure where black squares denote data present and white
%denotes not data present
set_figure_size(1);
pcolor(X,Y,freq); hold on; c = colormap('gray'); colormap(gca,flipud(c));
axis([1,ns+1,min(log10(1./T)),max(log10(1./T))+1])

for i=1:ns
 h = text(i,max(log10(1./T))+0.5,site{i});
 set(h,'Rotation',-90, 'fontsize',7, 'fontweight','bold');
end

ylabel('log10(Frequency)')
xlabel('Station Number')
title(['Frequency Analysis: ',num2str(length(uniqueT)),' distinct frequencies found. Max Freq = ',num2str(max(1./T)),'. Min Freq = ',num2str(min(1./T))])
set(gca,'Layer','top')
print_figure('.', 'Frequency_Analysis')

%Plot histogram of all frequencies present
%Note: Histograms for raw and interpolated data should look roughly the
%same
figure
histogram(log10(T_all))
title('Histogram of Periods Present');
xlabel('log10(Period)')
ylabel('Count')


   
  
  