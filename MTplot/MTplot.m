function MTplot
%
% Function which loads MT data and allows you to view apparent resistivity,
% phase, polar diagrams, phase tensors, tipper, induction vectors,
% dimensionless parameters, Bostick mapping and D+ fitting.
%
%
% Usage: MTplot
%   No inputs required
%
% There is the option to load a single EDI with no interpolation (e.g. load
% raw data), or to load multiple EDIs (all in current directory) and
% interpolate or bin the EDIs to a common frequency set. After doing the
% interpolation, the data is saved in a matfile which contains the "d" MT
% data structure. This matfile can then be loaded (rather than doing the
% interpolation or binning every time). There is also a "legacy" matfile
% from previous versions of MTplot_vxx which can also be loaded but this has
% not been thoroughly de-bugged.
%
%
% (C) 2020 Unsworth Research Group (University of Alberta, Edmonton, Canada)

clear all; close all;


curdir = pwd;

iload = menu('Load:','Single EDI',' Interpolate All EDI Files in Current Folder', 'Mat File (d structure)','Mat File (legacy)','ModEM Format File','Exit');

if iload == 1 %Option 1: Load a single EDI file (no interpolation)
    
    pause(0.1)
    [edifile,edipath] = uigetfile({'*.edi'},'Pick EDI file');
    
    if edifile==0
        return
    end
    
    cd(edipath);
    d = load_data_edi(edifile);
    cd(curdir);
    
elseif iload == 2 %Option 2: Load all edi files in folder
    
    d=interpolate_edi;

elseif iload == 3 %Option 3: Load matfile containing "d" data structure
    
    pause(0.1)
    [matfile,matpath] = uigetfile({'*.mat'},'Pick Mat File');
    
    if matfile==0
        return
    end
    
    cd(matpath);
    load(matfile);
    cd(curdir);
    
elseif iload == 4 %Option 4 = load matfile (LEGACY MATFILES)
    pause(0.1)
    [matfile,matpath] = uigetfile({'*.mat'},'Pick Mat File (Legacy)');
    
    if matfile == 0
        return
    end
    
    cd(matpath);
    mat = load(matfile);
    cd(curdir);
    
    %Need to rename variables into d data structure standard format
    d.loc(:,1) = mat.coord(:,2);
    d.loc(:,2) = mat.coord(:,1);
    d.loc(:,3) = mat.coord(:,3);
    d.site = mat.site;
    [d.T, index] = sort(mat.T);
    d.f = 1./d.T;
    d.nf = length(d.T);     d.nr = 4;       d.ns = length(d.loc(:,1));
    d.responses = {'ZXX','ZXY','ZYX','ZYY'};
    
    mu = 4*pi*10^-7;
    d.Z = (mu*1000)*permute(mat.ZZZ,[2 1 3]);
    d.Z = d.Z(index,:,:);
    d.Zerr = (mu*1000)*permute(mat.EEE,[2,1,3]);
    d.Zerr = d.Zerr(index,:,:);
    d.tip = permute(mat.ZZZ_HZ,[2 1 3]);
    d.tip = d.tip(index,:,:);
    d.tiperr = permute(mat.EEE_HZ,[2 1 3]);
    d.tiperr = d.tiperr(index,:,:);
    
    [d.rho, d.pha, d.rhoerr, d.phaerr] = calc_rho_pha(d.Z,d.Zerr,d.T);
    
    for is = 1:length(d.site)
        d.zrot(:,is) = mat.station(is).zrot*ones(d.nf,1);
        d.trot(:,is) = mat.station(is).trot*ones(d.nf,1);
    end
    
elseif iload == 5 %Load ModEM format file
    
    pause(0.1)
    [modemfile,modempath] = uigetfile({'*.dat;*.data', 'Data Files (*.dat,*.data)'; '*.*', 'All Files (*.*)'},'Pick ModEM Data File');
    
    if modemfile==0; return; end
    cd(modempath);
    d = load_data_modem(modemfile);
    cd(curdir);
    
else
    return
    
end

if nanmean(log10(d.rho(:)))>6 %Check if apparent resistivities are huge
    disp('WARNING: Your apparent resistivity values are very large. Your Z impedance values may not be in SI units.')
end

% Set up parameters for m_map and geoboundaries
u = user_defaults;
[d] = set_map_projection(d);
[L] = load_geoboundary_file_list;

%==========================================================================
% Start main plotting loop
irun = 1;
while irun == 1
    
    ichoice = menu('','Plot SRTM Station Map','Plot Resistivity and Phase','Plot Polar Diagrams','Plot Phase Tensor', ...
        'Plot Tipper, TFs, IVs','Plot Dimensionless Parameters','Plot Bostick Mapping (not working yet)',...
        'Edit Data','Reload user defaults','Quit');
    
    close all
    
    if ichoice == 1 %PLOT SRTM MAP*****************************************
        % Plot map of station locations with SRTM data
        plot_stations_map_srtm(d);
    
    elseif ichoice == 2 %APPARENT RESISTIVITY AND PHASE************************
        
        is = 1; rho_run = 1;
           
        while rho_run == 1
        rho_pha_menu = menu('Rho Pha','Plot All Stations At Once','Plot Individual Station','Pick Station From Map','Interpolated On Map','Pseudo-Section Profile','Return');

        if rho_pha_menu == 1 %PLOT ALL STATIONS AT ONCE----------------

            % Plot simple figure to check data
            [~] = set_figure_size(1);
            %Not sure if its really useful to plot rotations, probably more
            %useful to plot phase
%             subplot(221); plot(1:d.ns,d.zrot(1,:),'ob'); title('Rotations')
%             xlabel('Station number'); ylabel('Azimuth Rotation (deg)')

            for idxs=1:d.ns
                subplot(221); title('Apparent Resistivity')
                loglog(d.T,d.rho(:,2,idxs),'r'); hold on; loglog(d.T,d.rho(:,3,idxs),'b');
                
                subplot(223); title('Phase')
                semilogx(d.T,d.pha(:,2,idxs),'r'); hold on; semilogx(d.T,d.pha(:,3,idxs)+180,'b');
                xlabel('Period (s)');

                subplot(224); title('Tipper')
                semilogx(d.T,abs(d.tip(:,1,idxs)),'r'); hold on; semilogx(d.T,abs(d.tip(:,2,idxs)),'b');
                xlabel('Period (s)');
            end
            

            subplot(222);
            plot(1:d.ns,d.loc(:,3),'ob');
            title('Station Elevations')
            ylabel('Elevation (m)')
            xlabel('Station Number')

        elseif rho_pha_menu == 2 %PLOT STATIONS ONE-BY-ONE-------------

            dr = d; rotation = 0; dplus_flag = 0;
            
            while 1

                set_figure_size(1);
                plot_rho_pha(dr,is)
                mh = subplot(2,3,[3 6]);
                plot_topo(dr,1);
                plot_geoboundaries_geoshow(L)
                geoshow(d.loc(:,1),d.loc(:,2),'DisplayType','point','Marker','.','MarkerEdgeColor','k'); hold on
                geoshow(d.loc(is,1),d.loc(is,2),'DisplayType','point','Marker','.','MarkerEdgeColor','r','MarkerSize',15)
                axis(d.lim)
                xlabel('Longitude'); ylabel('Latitude')

                subplot(2,3,1)
                title(['Rotation = ',num2str(rotation),char(176)])
                
                if dplus_flag
                    [dpxy] = dplus(dr.Z(:,2,is),dr.Zerr(:,2,is),d.T);
                    [dpyx] = dplus(dr.Z(:,3,is)*exp(1i*pi),dr.Zerr(:,3,is),d.T);
                    
                    [~, ~, rhoerrx, phaerrx] = calc_rho_pha(dpxy.Z,dpxy.Zerr,dpxy.T);
                    [~, ~, rhoerry, phaerry] = calc_rho_pha(dpyx.Z,dpyx.Zerr,dpyx.T);
                    
                    nfig = length(findobj(gcf,'type','axes'));

                    if nfig == 5
                        subplot(2,3,2)
                    else
                        subplot(2,2,1)
                    end

                    loglog(dpxy.T,dpxy.rho,'-r','LineWidth',2);
                    loglog(dpyx.T,dpyx.rho,'-b','LineWidth',2);
                    
                    idx = ~isnan(dr.rho(:,2,is));
                    idy = ~isnan(dr.rho(:,3,is));
                    logerrorbar(dpxy.T,dr.rho(idx,2,is),rhoerrx,'or','-r');
                    logerrorbar(dpyx.T,dr.rho(idy,3,is),rhoerry,'sb','-b');
                    manual_legend('XY','or','YX','sb');

                    if nfig == 5
                        subplot(2,3,5)
                    else
                        subplot(2,2,3)
                    end

                    semilogx(dpxy.T,dpxy.pha,'-r','LineWidth',2);
                    semilogx(dpyx.T,dpyx.pha,'-b','LineWidth',2);
                    
                    errorbar(dpxy.T,dr.pha(idx,2,is),phaerrx,'or');
                    errorbar(dpyx.T,dr.pha(idy,3,is)+180,phaerry,'sb')

                    
                    title(mh,['TE R.M.S. = ',num2str(dpxy.rms),'. TM R.M.S. = ',num2str(dpyx.rms)])
                    
                end

                print_figure('rho_pha',[d.site{is},'_rot_',num2str(rotation)]); %Save figure

                next_menu = menu('','Next Station','Previous Station','Rotate +5 degrees','Rotate -5 degrees','D+ Fit On','D+ Fit Off','Return');

                if next_menu == 1

                    is = is+1;
                    if is > d.ns
                        is = d.ns;
                    end
                elseif next_menu == 2

                    is = is-1;
                    if is<1
                        is = 1;
                    end
                    
                elseif next_menu == 3 %Rotate 5 degrees clockwise
                    
                    rotation = rotation + 5;
                    dr = rotate_d(dr,5);
                    
                elseif next_menu == 4 %Rotate 5 degrees counter-clockwise
                    
                    rotation = rotation - 5;
                    dr = rotate_d(dr,-5);
                    
                elseif next_menu == 5 %D+
                    
                    dplus_flag = 1;
                    
                elseif next_menu == 6
                    
                    dplus_flag = 0;
                    
                else
                    break
                end

            end

        elseif rho_pha_menu == 3 %SELECT STATION FROM MAP--------------

            while 1
                figure(2)
                plot_topo(d,1);
                plot_geoboundaries_geoshow(L)
                geoshow(d.loc(:,1),d.loc(:,2),'DisplayType','point','Marker','.','MarkerEdgeColor','k'); hold on
                axis(d.lim)
                xlabel('Longitude'); ylabel('Latitude')
                [long,lat] = ginput(1);

                [~,is] = min(distance(long,lat,d.loc(:,2),d.loc(:,1)));
                geoshow(d.loc(is,1),d.loc(is,2),'DisplayType','point','Marker','.','MarkerEdgeColor','r','MarkerSize',15)

                set_figure_size(1);
                plot_rho_pha(d,is)
                subplot(2,3,[3 6])
                plot_topo(d,1);
                plot_geoboundaries_geoshow(L)
                geoshow(d.loc(:,1),d.loc(:,2),'DisplayType','point','Marker','.','MarkerEdgeColor','k'); hold on
                geoshow(d.loc(is,1),d.loc(is,2),'DisplayType','point','Marker','.','MarkerEdgeColor','r','MarkerSize',15)
                axis(d.lim)
                xlabel('Longitude'); ylabel('Latitude')

                %print_figure('rho_pha',[d.site{is},'_rot_',num2str(rotation)]); %Save figure

                next_menu = menu('','Select Another','Return');

                if next_menu ~= 1
                    break
                end
            end


        elseif rho_pha_menu == 4 %INTERPOLATE ON MAP-------------------
            
            if d.ns == 1
                disp('Only one station loaded. Not possible to interpolate')
            else

                det_menu = menu('','Plot all components','Plot determinant rho and pha','Return');
                %Choose period range
                per_range={[num2str(1),':',num2str(1),':',num2str(d.nf)]};
                answer=char(inputdlg({'Enter period range (ex. : 1:1:13)'},'Plot Periods',1,per_range));
                %The periods to plot (ptp = periods to plot)
                ptp=str2num(answer);  %#ok<*ST2NM>

                if isempty(ptp)
                    ptp = 1;
                    disp('Invalid Entry: Plotting Period 001 as Default');
                end

                if max(ptp) > d.nf || min(ptp) < 1
                    ptp = 1;
                    disp('Invalid Entry: Plotting Period 1 as Default')
                end

                for ip = ptp
                    if det_menu == 1
                        plot_rho_pha_map(d,ip)
                        annotation('textbox', [0.3 0.9 1 0.08], ...
                            'String', ['Rotation = ',num2str(d.zrot(1))], ...
                            'EdgeColor', 'none', ...
                            'HorizontalAlignment', 'center','FontSize',12,'FontWeight','bold')
                        print_figure('rho_pha_map',['rho_pha_map_',num2str(ip,'%03.0f'),'_',num2str(d.T(ip)),'_rot_',num2str(d.zrot(1))]); %Save figure
                    elseif det_menu == 2
                        plot_rho_pha_map_determinant(d,ip)
                        annotation('textbox', [0.3 0.9 1 0.08], ...
                            'String', ['Rotation = ',num2str(d.zrot(1))], ...
                            'EdgeColor', 'none', ...
                            'HorizontalAlignment', 'center','FontSize',12,'FontWeight','bold')
                        print_figure('rho_pha_map_det',['rho_pha_map_det_',num2str(ip,'%03.0f'),'_',num2str(d.T(ip)),'_rot_',num2str(d.zrot(1))]); %Save figure
                    end

                end
            end

        elseif rho_pha_menu == 5 %PSEUDO-SECTIONS----------------------
            
            if d.ns == 1
                disp('Only one station loaded. Not possible to interpolate')
            else
                plot_rho_pha_pseudo(d);
            end

        else
            rho_run = 0;

        end
        
        end
        
    elseif ichoice == 3 %***************POLAR DIAGRAMS*********************
        
        
        
        while 1

            polar_menu = menu('Polar Diagrams','By Station','On Map','Distort Polar Diagrams','Return');

            if polar_menu == 1 %By station
                is = 1; %Initialize station counter
                
                while 1
                    s = d;
                    s.Z = d.Z(:,:,is);   s.Zerr = d.Zerr(:,:,is);  
                    s.tip = d.tip(:,:,is);    s.tiperr = d.tiperr(:,:,is); 
                    s.ns = 1;
                    s.rho = d.rho(:,:,is); s.rhoerr = d.rhoerr(:,:,is);
                    s.pha = d.pha(:,:,is); s.phaerr = d.phaerr(:,:,is);
                    s.name = [d.name,'_Site_',d.site{is}];
                    s.zrot = d.zrot(:,is); s.trot = d.trot(:,is);
                    s.site{1} = d.site{is}; s.loc = d.loc(is,:);
                    [polar] = calc_polar(s);
                    plot_polar(s,polar)
                    
                    print_figure('polar',['polar_station_',num2str(is,'%03.0f'),'_',d.site{is}]); %Save figure

                    next_menu = menu('','Next Station','Previous Station','Return');

                    if next_menu == 1
                        is = is+1;
                        if is>d.ns
                            is = d.ns;
                        end
                    elseif next_menu == 2
                        is = is-1;
                        if is<1
                            is = 1;
                        end
                    else
                        break
                    end
                
                end
               
            elseif polar_menu == 2 %On Map
                
                disp('Computing Polar Diagrams for all stations and all frequencies (Takes about 5 - 10 seconds)')
                [polar] = calc_polar(d);
                plot_polar_map(d,polar)

            elseif polar_menu == 3 %Distort polar diagrams (does not work yet)
                
                plot_polar_distorted(d)
            else
                break

            end
        end
        
    elseif ichoice == 4 %***************PHASE TENSOR***********************
        
        while 1
            phase_tensor_menu = menu('Phase Tensor','By Station','On Map','On Map With IVs','Ellipse Profile Pseudo-Section','Pseudo-Section Phimin Phimax','Beta Skew / Eccentricity /Strike Rose Diagrams','Return');

            if phase_tensor_menu == 1 %Phase Tensor By station

                plot_phase_tensor(d);

            elseif phase_tensor_menu == 2 %Phase Tensor on Map

                plot_phase_tensor_map(d,0);

            elseif phase_tensor_menu == 3 %Phase Tensor on Map with induction vectors
                
                plot_phase_tensor_map(d,1);

            elseif phase_tensor_menu == 4 %Phase Tensor Ellipse Profile
                
                close all
                plot_phase_tensor_ellipse_pseudo(d);

            elseif phase_tensor_menu == 5 %Phimin Phimax Pseudo Section
                
                close all
                plot_phase_tensor_pseudo(d);
                
            elseif phase_tensor_menu == 6 %Phase Tensor Parameters
                
                plot_phase_tensor_beta_strike(d)
                
            else
                break

            end
        
        end
            
        
    elseif ichoice == 5 %*************TIPPER AND INDUCTION VECTORS*********
        
        hz_menu = menu('Tipper, TFs, IVs','TFs and IVs by Station','IV on Map','Tipper Pseudo Section','IV Interpolated Angle','Return');
        
        if hz_menu == 1 %Transfer Functions and IVs by Stations
            
            is = 1;
            while 1
                
                set_figure_size(1);
                plot_tipper(d,is);
                plot_induction_vector_site(d,is);
                
                next_menu = menu('','Next Station','Previous Station','Quit');
                
                if next_menu == 1
                    is = is+1;
                    if is>d.ns
                        is = d.ns;
                    end
                elseif next_menu == 2
                    is = is-1;
                    if is<1
                        is = 1;
                    end
                else
                    break
                end
            end
            
            
        elseif hz_menu == 2 %Induction Vectors on Map
            
            for ip = 1:u.nskip:d.nf
                set_figure_size(1);
                plot_induction_vector_map(d,ip,'k');
                title(['Induction Vectors for Period: ',num2str(d.T(ip)),' s'])
                print_figure('iv_map',['iv_map_',num2str(ip,'%03.0f'),'_',num2str(d.T(ip))]); %Save figure
                pause(0.1)
            end
            
        elseif hz_menu == 3 %Real and Imag Induction Vector Pseudo Section
            
            plot_tipper_pseudo(d);
            
        elseif hz_menu == 4 % Induction Vector Angle Interpolated on Map
            
            plot_tipper_angle_map(d)
            
        end
        
    elseif ichoice == 6 %***********DIMENSIONLESS PARAMETERS***************
         
        dim_menu = menu('','By Station','Skew Pseudo Section');
        
        if dim_menu == 1
            plot_dim_parameters(d);
        elseif dim_menu == 2
            plot_skew_pseudo(d);
        end
        
    elseif ichoice == 7 %*************BOSTICK******************************
        
        bostick_menu = menu('Bostick','By Station','On Map','Profile');
        
        if bostick_menu == 1
            
            
        elseif bostick_menu == 2
            
            
        elseif bostick_menu == 3
            
        end
  
    elseif ichoice == 8 %********************EDIT DATA********************
        
        edit_menu = menu('','Edit Apparent Resisitivity','Edit Tipper','Return');
        
        if edit_menu==1
            
            d = edit_data(d);
            
        elseif edit_menu==2
            
            d = edit_data_tipper(d);
            
        end
            
    elseif ichoice == 9 %****************RELOAD USER DEFAULTS**************
        
        clear 'user_defaults.m' % this is needed to load changes saved in editor while S3D is running
        u = user_defaults;
        
        u_path = which('user_defaults.m');       
        disp(['Reloaded user_defaults from ',u_path])
                                   
    else
        break

    end
    
end
    
