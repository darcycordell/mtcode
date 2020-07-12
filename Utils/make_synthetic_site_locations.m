function d = make_synthetic_site_locations
% Function which allows "empty" data structure to be created with
% station locations for use in synthetic studies. These station locations 
% can then be used to create a model mesh in M3D.
%
% Option to load station locations from text file, matfile (d structure),
% ModEM data file, manually create a uniform grid, manually create a
% staggered grid, or randomly distribute stations in a rectangular area
%
% Outputs a matfile containing "d" data strucutre which can be loaded into
% M3D.
%
% The "empty" data structure has data with all 1's. Elevation is zero
% unless loaded stations contain elevation data.
%
% Note there is no way to add topography to synthetic data yet unless the
% data are loaded from a textfile or matfile which already contains
% topography. So if you want to add topography you will have to manually
% set station elevations in the data structure.
%

clear all
curdir = pwd;

site_menu = menu('Get Site Locations From:','Text File (lat,long,z)','Matfile "d" Structure','ModEM Data File','Manual Uniform Array (Default)','Manual Staggered Array','Random Site Locations');

if site_menu == 1 %TEXT FILE
    [d] = make_ones_data;
    
    [site_file,site_path] = uigetfile({'*.txt'},'pick  file containing site locations (lat,long,z)');
    cd(site_path)
    sites = load(site_file);
    cd(curdir)
    sze = size(sites);
    d.loc = [sites(:,1) sites(:,2)];
    
    if sze(2) == 2
        d.z = lat*0;
    else
        d.z = sites(:,3);
    end
    
    d.loc(:,3) = d.z;
    
    d.origin = [(max(d.loc(:,1))-min(d.loc(:,1)))/2+min(d.loc(:,1)),(max(d.loc(:,2))-min(d.loc(:,2)))/2+min(d.loc(:,2)),0];
    [d.y,d.x] = geo2utm(d.loc(:,2),d.loc(:,1),d.origin(2),d.origin(1));
    d.y = d.y-500000;
        
elseif site_menu == 2 %MATFILE
    [d] = make_ones_data;
    
    [mat_file,mat_path] = uigetfile({'*.mat'},'pick  matfile containing site locations (d structure)');
    cd(mat_path)
    dmat = load(mat_file);
    dmat = dmat.d;
    cd(curdir)
    
    d.loc = dmat.loc;
    
    d.origin = [(max(d.loc(:,1))-min(d.loc(:,1)))/2+min(d.loc(:,1)),(max(d.loc(:,2))-min(d.loc(:,2)))/2+min(d.loc(:,2)),0];
    [d.y,d.x] = geo2utm(d.loc(:,2),d.loc(:,1),d.origin(2),d.origin(1));
    d.y = d.y-500000;
    d.z = d.loc(:,3);
    
elseif site_menu == 3 %MODEM DATA FILE
    [d] = make_ones_data;
    
    [mod_file,mod_path] = uigetfile({'*.data';'*.dat'},'pick  ModEM data file containing site locations');
    cd(mod_path)
    dmat = load_data_modem(mod_file);
    cd(curdir)
    
    d.loc = dmat.loc;
    d.origin = dmat.origin;
    d.y = dmat.y;
    d.x = dmat.x;
    d.z = dmat.z;
    
elseif site_menu == 6 %RANDOM ARRAY
    prompt={'Minimum Station Spacing (km)','Total Number of Stations','Array Dimensions [minx maxx miny maxy] (km)'};
    dlg_title='Random Array of Stations';
    def={'4','60','[-40 40 -40 40]'};
    num_lines=1;
    dinp = inputdlg(prompt,dlg_title,num_lines,def);
    
    if isempty(dinp)
        return
    end
    
    ds = str2double(dinp{1})*1000;
    ns = str2double(dinp{2});
    g = str2num(dinp{3});
    minx = g(1)*1000;
    maxx = g(2)*1000;
    miny = g(3)*1000;
    maxy = g(4)*1000;
    A = (maxx-minx)*(maxy-miny)/10^6;
    
    tot_stn = 0; count = 0;
    x = rand(ns,1)*2*maxx+minx;
    y = rand(ns,1)*2*maxy+miny;
    while tot_stn<ns
        
        tot_stn = length(x);
        dist = zeros(tot_stn); idx = zeros(tot_stn); 
        for i = 1:length(x)
            dist(i,:) = sqrt((x(i)-x).^2+(y(i)-y).^2);
            idx(i,:) = dist(i,:)<ds;
            idx(i,i) = 0;
        end

        for i = 1:length(x)
            if any(idx(:,i))
                x(i) = NaN;
                y(i) = NaN;
            end
        end

        x(isnan(x))=[];
        y(isnan(y))=[];
        tot_stn = length(x);
        
        x = [x; rand(ns-tot_stn,1)*2*maxx+minx];
        y = [y; rand(ns-tot_stn,1)*2*maxy+miny];
        
        count = count+1;
        
        if count>20
            count = 0;
            ds = ds/1.1;
        end

    end
    
    if ds~=str2double(dinp{1})*1000
        disp(['***It is not possible to have ',num2str(ns),' sites randomly distributed in a ',num2str(A),' sq km area with a station spacing of ',dinp{1},' km.'])
        disp(['The minimum station spacing has been optimized to be ',num2str(ds/1000),' km***'])
    end
    
    [d] = make_ones_data;
    
    [long,lat] = utm2geo(y,x,0,0);

    d.loc = [lat long x.*0];
    d.x = x;
    d.y = y;
    d.z = x.*0;


else %MANUAL UNIFORM GRID OR MANUAL STAGGERED GRID
    prompt={'NS Station Spacing (km)','EW Station Spacing (km)','Station Array Dimensions [minx maxx miny maxy] (km)'};
    dlg_title='Grid of Stations';
    def={'4','4','[-40 40 -40 40]'};
    num_lines=1;
    dinp = inputdlg(prompt,dlg_title,num_lines,def);
    
    if isempty(dinp)
        return
    end
    
    dx = str2double(dinp{1})*1000;
    dy = str2double(dinp{2})*1000;
    g = str2num(dinp{3});
    minx = g(1)*1000;
    maxx = g(2)*1000;
    miny = g(3)*1000;
    maxy = g(4)*1000;
    
    x = [minx:dx:maxx];
    y = [miny:dy:maxy];

    [d] = make_ones_data;
    [X,Y]=meshgrid(x,y);
    
    if site_menu == 5 %STAGGER THE GRID
        X = X(1:2:end);
        Y = Y(1:2:end);
    end
            
    x = X(:);
    y = Y(:);
    z = x.*0;

    [long,lat] = utm2geo(y,x,0,0);

    d.loc = [lat long z];
    d.x = x;
    d.y = y;
    d.z = z;
    
end

freq_menu = menu('Get Frequency List From:','Frequency Text File','Manually Enter (Default)');

if freq_menu == 1 %LOAD FREQUENCY TEXT FILE
    [freq_file,freq_path] = uigetfile({'*.txt'},'pick  file containing frequencies');
    cd(freq_path)
    f = load([freq_path freq_file],'-ascii');
    cd(curdir)
    
else %MANUALLY CREATE FREQUENCY LIST
    prompt={'Minimum Frequency (Hz)','Maximum Frequency (Hz)','Number of Frequencies'};
    dlg_title='Frequency Setup';
    def={'0.001','1000','40'};
    num_lines=1;
    dinp = inputdlg(prompt,dlg_title,num_lines,def);
    
    if isempty(dinp)
        return
    end
    
    minf = str2double(dinp{1});
    maxf = str2double(dinp{2});
    nf = str2double(dinp{3});
    
    f = 10.^linspace(log10(minf),log10(maxf),nf).';

end

d.f = flipud(f);
d.T = 1./d.f;
d.nf = numel(d.f);
d.responses = {'ZXX' 'ZXY' 'ZYX' 'ZYY'};
d.ns = length(d.x);

d.Z = repmat(d.Z,[d.nf 1 d.ns]);
d.Zerr = repmat(d.Zerr,[d.nf 1 d.ns]);
d.rho = repmat(d.rho,[d.nf 1 d.ns]);
d.rhoerr = repmat(d.rhoerr,[d.nf 1 d.ns]);
d.pha = repmat(d.pha,[d.nf 1 d.ns]);
d.phaerr = repmat(d.phaerr,[d.nf 1 d.ns]);
d.tip = repmat(d.tip,[d.nf 1 d.ns]);
d.tiperr = repmat(d.tiperr,[d.nf 1 d.ns]);

d.zrot = zeros(d.nf, d.ns);
d.trot = zeros(d.nf, d.ns);

d.site = cell(d.ns,1);
for i = 1:d.ns
    d.site{i} = ['Site_',num2str(i,'%03.0f')];
end

d = set_map_projection(d);

d.origin(1) = min(d.loc(:,1)) + (max(d.loc(:,1))-min(d.loc(:,1)))/2;
d.origin(2) = min(d.loc(:,2)) + (max(d.loc(:,2))-min(d.loc(:,2)))/2;
d.origin(3) = min(d.z);

prompt={'File Name'};
dlg_title='Save File';
def={'Synthetic_Site_Locations'};
num_lines=1;
dinp = inputdlg(prompt,dlg_title,num_lines,def);

d.name = dinp{1};

save(d.name,'d')

figure
subplot(1,2,1)
plot(d.y,d.x,'ro')
xlabel('y (m)'); ylabel('x (m)'); title('Synthetic Station Locations in model coordinates')
axis equal
subplot(1,2,2)
geoshow(d.loc(:,1),d.loc(:,2),'displaytype','point')
xlabel('Longitude'); ylabel('Latitude'); title('Synthetic Station Locations in geographic coordinates')