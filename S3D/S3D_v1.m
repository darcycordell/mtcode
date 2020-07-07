function S3D_v1(varargin)
% Function to plot 3D results from ModEM or WSINV Inversions.
%
% In order to run this code for ModEM, you will need three input files in your **current directory**:
%
% 1) The initial starting model file (usually a halfspace)
% 2) The *.rho 3D inversion model
% 3) The observed data file
% 4) The *.dat 3D inversion response file which is loaded automatically
% given that it has the same name as the *.rho file.
%
% These files can be specified in a parameter text file with three lines corresponding
% to the filenames (in the order above). Example: S3.par
%
% halfspace.model
% ModEM_NLCG_itr.rho
% observed_data.data
%
% In order to run this code for WSINV, you will need the inversion model
% file, the original mat file produced from load_edis_vxx, the inversion
% startup file and an option text file which includes the geological or
% geographic boundaries to be plotted on maps.
%
% wsinv_inversion_model.10
% data.mat
% startup
%
%
% The code loads the observed data, predicted data, starting model, and
% inversion model into the standard structure formats (dobs, dpred, m, and
% m0)
%
% Potential Bugs or Things to Add:-----------------------------------------
%
% -The code has not been thoroughly tested with WSINV3D
% -The code still uses m_map throughout. It might be useful to transition
% towards the MATLAB geo package.
%
% (C) 2020 Unsworth Research Group (University of Alberta, Edmonton, Canada)


close all

%------------------User input structures-----------------------------------
% This is a separate matlab function which contains a variety of useful
% user inputs such as axis limits, geoboundary files, colormaps etc.
% This function should be kept on the matlab path and each person has a
% different version.
u = user_defaults;
u_path = which('user_defaults.m');
disp(['Loaded user_defaults from ',u_path])

%--------------------------------------------------------------------------

%
%-----------------------BEGIN LOADING DATA---------------------------------
inversion_algorithm_menu = menu('','ModEM isotropic','ModEM anisotropic','WSINV3D','Quit');

if inversion_algorithm_menu == 1 || inversion_algorithm_menu == 2 %ModEM inversion

    %Get Filenames
    if isempty(varargin)==0 %If a parfile is supplied, then read it.

        parfile = varargin{1};

        fid = fopen(parfile);
        if fid==-1
            error([parfile,' parameter file not found'])
        else
            reading = textscan(fid,'%s');
            reading = reading{1};
        end
        fclose all;

        m0_file = reading{1}; %Input model filename
        m_file = reading{2}; %Inversion model filename
        dobs_file = reading{3}; %Observed data filename

    else %If no parfile is supplied, then get user input to choose files

        [m0_file]=uigetfile({'*.model'},'Choose Initial ModEM model file'); if m0_file == 0; return; end
        [m_file]=uigetfile({'*.rho'},'Choose an output model file from the ModEM inversion'); if m_file == 0; return; end
        [dobs_file]=uigetfile({'*.data'},'Choose Initial ModEM data file'); if dobs_file == 0; return; end

    end

    %Get other file names based off the inversion model file name
    if inversion_algorithm_menu == 1 % ModEM isotropic
        basefile=m_file(1:(length(m_file)-4));% base file name used for ModEM iteration file names  
    else % ModEM anisotropic, need to change how basefile name is read
        basefile=m_file(1:(length(m_file)-6));% base file name used for ModEM iteration file names  
        disp('anisotropic model loaded') % edit basefile name -4 to -6
    end
    
    logfile = [basefile(1:end-4) '.log']; %ModEM log file
    dpred_file = [basefile '.dat']; %ModEM predicted (or modelled) inversion response data file

    %Load Data
    [dobs] = load_data_modem(dobs_file);
    [dpred] = load_data_modem(dpred_file);

    %Load Model
    [m] = load_model_modem(m_file);
    [m0] = load_model_modem(m0_file);

    %Load logfile
    [invlog] = load_logfile_modem(logfile);

elseif inversion_algorithm_menu == 3 %WSINV3D inversion
    %Get Filenames
    if isempty(varargin)==0 %If a parfile is supplied, then read it.

        parfile = varargin{1};

        fid = fopen(parfile);
        if fid==-1
            error([parfile,' parameter file not found'])
        else
            reading = textscan(fid,'%s');
            reading = reading{1};
        end
        fclose all;

        m_file = reading{1}; %Inversion model filename
        matfile = reading{2}; %Matfile of observed data
        startupfile = reading{3}; %Startup filename

    else %If no parfile is supplied, then get user input to choose files

        [m_file]=uigetfile({'*'},'Choose an output model file from the WSINV inversion'); if m_file == 0; return; end
        [matfile]=uigetfile({'*.mat'},'Choose MAT data file'); if matfile == 0; return; end
        [startupfile]=uigetfile({'*'},'Choose the inversion startup file'); if startupfile == 0; return; end

    end
    
    %Read the startup file to get other filenames
    [invtype, dobs_file, ~, ~, logfile, m0_file] = load_startup_wsinv(startupfile);
%     fid = fopen(startupfile,'r');
%     tmp=0;
%     while isempty(tmp)~=1
%         tmp = fscanf(fid,'%s',1); % reading
%         if strcmp(tmp,'INVERSION_TYPE')
%             invtype = str2double(fscanf(fid,'%s',1));
%         elseif strcmp(tmp,'DATA_FILE')
%             dobs_file = fscanf(fid,'%s',1);
%         elseif strcmp(tmp,'OUTPUT_FILE')
%             logfile = [fscanf(fid,'%s',1),'.log'];
%         elseif strcmp(tmp,'INITIAL_MODEL_FILE')
%             m0_file = fscanf(fid,'%s',1);
%         end
%     end

    %Get other file names based off the inversion model file name
    niter = m_file(end-1:end);
    basefile=m_file(1:(length(m_file)-9));% base file name used for WSINV iteration file names  
    dpred_file = [basefile '_resp.',niter]; %WSINV predicted (or modelled) inversion response data file

    %Load Data
    [dobs] = load_data_wsinv(dobs_file,matfile,invtype);
    [dpred] = load_data_wsinv(dpred_file,matfile,invtype);

    %Load Model
    [m] = load_model_wsinv(m_file);
    [m0] = load_model_wsinv(m0_file);

    %Load logfile
    [invlog] = load_logfile_wsinv(logfile);

else
    return
end

%Check if geofile exists and if not, set geofile to 'none'
fid = fopen(u.geofile);
if fid== -1 %If no geo file is found on path, then set geofile to none and it will not be plotted
    disp('Geoboundary File cannot be found on your path or was not specified and so will not be plotted')
    u.geofile = 'none';
end

%-----------------DONE LOADING DATA----------------------------------------

%--------------------------------------------------------------------------

%Add to the data and model structures all the necessary conversions from
%model coordinates to latitude and longitude
disp('Converting Model Space to Lat-Long and Data Coordinates to UTM')
[m,dobs] = link_model_data(m,dobs); %Add lat-long to inversion model

disp('Converting Input Model Space to Lat-Long and Data Coordinates to UTM')
[m0,~] = link_model_data(m0,dobs); %Add lat-long to initial model

dpred.x = dobs.x; dpred.y = dobs.y;
m0.cy = m.cy; m0.cx = m.cx; m0.origin = m.origin; m0.x = m.x; m0.y = m.y;
m0.X = m.X; m0.Y = m.Y;

disp('Converting Predicted Data Coordinates to UTM')
[~,dpred] = link_model_data(m,dpred); %Add x-y to predicted data


%Option to plot either the input model or the inversion model
%"niter" variable is set accordingly and is used to specify which inversion
%iteration is being plotted.
if u.plot_input_model
    dobs.niter = 'input'; dpred.niter = 'input';
    m0.niter = 'input';
    model_to_plot = m0; %Model to plot is the input model
else
    %Get the inversion iteration number (different for WSINV vs. ModEM)
    if inversion_algorithm_menu == 1
        niter=basefile((length(basefile)-2):end); 
    elseif inversion_algorithm_menu == 2
        niter=['0',basefile((length(basefile)-1):end)];
    end
    
    dobs.niter = niter; dpred.niter = niter;
    m.niter = niter;
    
    model_to_plot = m; %Model to plot is the inversion model
end


[L] = load_geoboundary_file_list;

%------------------------BEGINNING MAIN MENU-------------------------------
menu1 = 1;

while menu1 == 1
    
    main_menu = menu('Main Menu','View Model','View Data/Response','Edit Model','Edit/Export Data','Reload user defaults','Quit');

    if main_menu == 1 %VIEW MODEL------------------------------------------------------------------------------
        menu2 = 1;

        while menu2 == 1
            
            model_menu = menu('','Vertical Cross-Section','Horizontal Slice (Lat/Lon)','Horizontal Slice (model X,Y)','Conductance','Rho-Depth Curves','Isosurface','3D View','Compare with Another Model','Export All Slices in Lat-Long','Plot Slices in 3D','Return');

            if model_menu==1 %Plot Vertical Cross-Sections
                close all
                while 1
                    section_menu = menu('','NS Sections','EW Section','Diagonal Section','Section Beneath Specified Stations','Export Current Cross Section in Lat-Long','Return');
                
                    if section_menu == 1 %Plot NS section----------------------
                        [x,y,z,rho] = plot_cross_section_3d(model_to_plot,'NS',dobs);     
                    elseif section_menu == 2 %Plot EW section------------------
                        [x,y,z,rho] = plot_cross_section_3d(model_to_plot,'EW',dobs);
                    elseif section_menu == 3 %Plot Diagonal Section------------
                        if strcmp(u.diagonal_section_mode,'line') || strcmp(u.diagonal_section_mode,'interp')
                            %This version takes line segment distances, no
                            %smoothing. Darcy's version.
                            [x,y,z,rho] = plot_diagonal_section(model_to_plot,dobs);
                        elseif strcmp(u.diagonal_section_mode,'smooth')
                            %This version averages across the line segment and
                            %results in a smoothed profile. Enci's version.
                            [x,y,z,rho] = plot_diagonal_section_smooth(model_to_plot,dobs);
                            
                        end

                    elseif section_menu == 4 %Export section
                        %Plots a fence diagram beneath selected stations
                        [x,y,z,rho]=plot_cross_section_beneath_points(m,dobs);

                    elseif section_menu == 5 %Export profile beneath stations
                        
                        g=groot;
                        if isempty(g.Children)
                            disp('In order to export a N-S, E-W, or diagonal section, you must first choose to plot one')
                        else
                            %Write out the section in latitude and
                            %longtiude.
                            write_section(x,y,z,rho,dobs);
                            
                            %Plot the section that you just output (for
                            %debugging purposes)
                            %plot_cross_section_lat_long
                        end

                    else
                        break

                    end
                
                end
                

            elseif model_menu == 2 %Plot Horizontal Slices in Lat/Lon-----------------
                prompt = {'Choose slice number(s) e.g.  1:25 or single number'};
                titles  = 'Choose slice';
                def = {['1:',num2str(model_to_plot.nz)]};
                slices = str2num(char(inputdlg(prompt,titles,1,def))); 

                %Loop over z-slices
                for i=1:length(slices)
                    close all
                    plot_slice_map(model_to_plot,slices(i),dobs);                    
                    title(['Depth = ',num2str(model_to_plot.cz(slices(i))/1000),' km b.s.l.']);
                    print_figure(['horizontal_slices_geo_',model_to_plot.niter],['slice_',num2str(slices(i),'%03.0f')]) 
                end
                
            elseif model_menu == 3 %Plot Horizontal Slices in model coordinates (X,Y)-----------------
                prompt = {'Choose slice number(s) e.g.  1:25 or single number'};
                titles  = 'Choose slice';
                def = {['1:',num2str(model_to_plot.nz)]};
                slices = str2num(char(inputdlg(prompt,titles,1,def))); 

                %Loop over z-slices
                for i=1:length(slices)
                    close all
                    plot_slice(model_to_plot,slices(i),dobs);                    
                    title(['Depth = ',num2str(model_to_plot.z(slices(i))/1000),' km b.s.l.']);
                    print_figure(['horizontal_slices_xy_',model_to_plot.niter],['slice_',num2str(slices(i),'%03.0f')]) 
                end

            elseif model_menu == 4 %Plot Conductance-----------------------
                
                plot_conductance(model_to_plot,dobs)
                print_figure(['horizontal_slices_',model_to_plot.niter],'conductance_map') 
                    
            elseif model_menu == 5 %Plot resistivity and depth curves------
                
                 rho_z = menu('','View Rho-z Curve Beneath Stations','Pick Model Cell to plot Rho-z Curve','Back');

                 if rho_z ==1
                    %View the resitivity-depth curve beneath stations only 
                    plot_rho_z_station(model_to_plot,dobs)
                    
                 elseif rho_z == 2
                    %View the resistivity-depth curve beneath any model
                    %cell (more general)
                    plot_rho_z_model(model_to_plot,dobs)
                    
                 end
                 
            elseif model_menu == 6 %Plot isosurface------------------------
                close all
                plot_isosurface(m,dobs)
                
            elseif model_menu == 7 %Plot 3D View---------------------------
                plot_model_3D(m,dobs)
                
            elseif model_menu == 8 %Compare Currently loaded Model to another model---------
                
                close all
                %This is currently only implemented for loading ModEM
                %models
                
                %Also, note that if one of your models has topography
                %while the other does not it is not possible to properly
                %compare the models.
                
                %Load a second model
                curdir = pwd;
                [filename,modpath]=uigetfile({'*.rho'},'Choose second ModEM model to compare'); if filename == 0; return; end
                cd(modpath);
                [model_to_compare] = load_model_modem(filename);
                
                [filename,datpath]=uigetfile({'*.dat'},'Choose datafile which corresponds to the model you chose'); if filename == 0; return; end
                cd(datpath);
                [data_to_compare] = load_data_modem(filename);
                
                cd(curdir);

                compare_models(model_to_plot,model_to_compare,dobs,data_to_compare);
                
            elseif model_menu == 9 %Export slices and cube in lat-long----------------------
                
                write_slices_cube_lat_long(model_to_plot);
                
            elseif model_menu == 10 % Plot slices in 3D--------------------------------------
                
                plot_slice_multi(m,dobs);
                
            else
                menu2 = 0;

            end

        end

    elseif main_menu == 2 %VIEW DATA------------------------------------------------------------------------------
        menu3 = 1;

        while menu3 == 1
            
            data_main_menu = menu('Data/Response Menu','Plot Rho / Phase Curves','Plot Impedance Curves','Plot Tipper Curves','Plot Rho / Phase Interpolated Map','Plot Induction Vector Map','Plot Rho / Phase Pseudo','Plot Tipper Pseudo','Misfit Statistics','Compare to another data response','K-S test with a third data response','Residuals cross-plot with a third data response','Return');

            if data_main_menu == 1 || data_main_menu == 2 || data_main_menu == 3 %If you want to plot data as curves------
                
                    s = detailed_statistics(dobs,dpred); % new 4/2020
                    plot_data(dobs,data_main_menu,dpred,s);
      
            elseif data_main_menu == 4 || data_main_menu == 5 %Plot Data in Map View-------------------------------
                    
                    %Choose period range
                    per_range={[num2str(1),':',num2str(1),':',num2str(dobs.nf)]};
                    answer=char(inputdlg({'Enter period range (ex. : 1:1:13)'},'Plot Periods',1,per_range));
                    %The periods to plot (ptp = periods to plot)
                    ptp=str2num(answer);  %#ok<*ST2NM>
                    
                    if isempty(ptp)
                        ptp = 1;
                        disp('Invalid Entry: Plotting Period 001 as Default');
                    end

                    if max(ptp) > dobs.nf || min(ptp) < 1
                        ptp = 1;
                        disp('Invalid Entry: Plotting Period 1 as Default')
                    end

                    %Loop over periods to plot
                    for ip = ptp

                        if data_main_menu == 4 %Plot interpolated rho and phase in map view
                            plot_rho_pha_map(dobs,ip);
                            print_figure(['compare_responses_map_',dpred.niter],['Observed','_',num2str(dobs.T(ip))]); %Save figure
                            plot_rho_pha_map(dpred,ip);
                            print_figure(['compare_responses_map_',dpred.niter],['Modelled','_',num2str(dpred.T(ip))]); %Save figure
                            plot_misfit_percent_diff_map(dobs,dpred,ip)
                        elseif data_main_menu == 5 %Plot induction vectors in map view
                            close all
                            if ~all(isnan(dobs.tip(:))) %If there is tipper data
                                plot_induction_vector_map(dobs,ip,'r')
                                plot_induction_vector_map(dpred,ip,'k')
                                plot_geoboundaries(L);
                                figure_name = 'induction_vector';
                                manual_legend('Observed Data','-r','Modelled Data','-k');
                                title(['Induction Vectors for T = ',num2str(dobs.T(ip))])
                                print_figure(['compare_responses_map_',dpred.niter],[figure_name,'_',num2str(dobs.T(ip))]); %Save figure
                            else
                                disp('No Tipper data in your observed data!')
                            end
                        end
                    end
                    
            elseif data_main_menu == 6 || data_main_menu == 7 % Plot Rho/Phase Pseudo section
                pseudo_menu = menu('Plot pseudo-section','Plot Observed Data','Plot Predicted Data','Return');
                
                if pseudo_menu == 1 % plot observed data
                    if data_main_menu == 6  % Plot Rho/Phase Pseudo section
                        plot_rho_pha_pseudo(dobs)
                    else % plot tipper pseudo
                        plot_tipper_pseudo(dobs)
                    end
                elseif pseudo_menu == 2 % plot predicted data
                    if data_main_menu == 6  % Plot Rho/Phase Pseudo section
                        plot_rho_pha_pseudo(dpred)
                    else
                        plot_tipper_pseudo(dpred)
                    end
                else
                    menu3 = 0;
                end

            elseif data_main_menu == 8 %Misfit Statistics----------------------------------------------

                misfit_menu = menu('','Residual Histograms','Misfit By Period','Misfit Map','Misfit Map (By Period)','Misfit Pseudo','Misfit By Inversion Iteration');

                s = detailed_statistics(dobs,dpred); % new 4/2020
                if misfit_menu == 1 %Residual histograms-------------------
                    close all                   
                    plot_misfit_residual_histograms(s)
                    print_figure(['misfit_stats_',dpred.niter],'Residuals'); %Save figure

                elseif misfit_menu == 2 %Misfit by period------------------
                    close all
                    plot_misfit_rms_by_period(dobs,dpred,s)
                    
                elseif misfit_menu == 3 %Misfit map view-------------------
                    close all;
                    plot_misfit_rms_map(dobs,dpred,s)    
                    plot_geoboundaries(L);     
                    print_figure(['misfit_stats_',dpred.niter],'rms_map_view'); %Save figure
                    
                elseif misfit_menu == 4 %Misfit map view (by period)
                    
                    close all
                    plot_misfit_rms_by_period_map(dobs,dpred,s)
                    
                elseif misfit_menu == 5 %Misfit pseudo-section  %% new 6/2020
                    close all
                    plot_misfit_residual_pseudo(dobs,dpred)
                    
                elseif misfit_menu == 6 %Misfit convergence----------------
                    plot_misfit_convergence(invlog);
                    print_figure(['misfit_stats_',dpred.niter],'inversion_convergence'); %Save figure
                end
                
            elseif data_main_menu == 9 % Compare to a third data response----------------------------------------------

                close all
                %This is currently only implemented for loading ModEM data               
                %Load a third data set
                curdir = pwd;               
                [filename,datpath]=uigetfile({'*.dat;*.data', 'ModEM data files (*.dat,*.data)'; '*.*', 'All Files (*.*)'},'Choose a ModEM data file to compare to'); if filename == 0; return; end
                cd(datpath);
                [data_to_compare] = load_data_modem(filename);            
                cd(curdir);

                compare_data(dobs,dpred,data_to_compare)        
                
            elseif data_main_menu == 10 % K-S test with a third data response----------------------------------------------

                close all
                %This is currently only implemented for loading ModEM data               
                %Load a third data set
                curdir = pwd;               
                [filename,datpath]=uigetfile({'*.dat;*.data', 'ModEM data files (*.dat,*.data)'; '*.*', 'All Files (*.*)'},'Choose a ModEM data file to compare to'); if filename == 0; return; end
                cd(datpath);
                [data_to_compare] = load_data_modem(filename);   
                cd(curdir);
                
                disp(['Using significance level of ',num2str(u.significance_level),' defined in user_defaults.'])
                [ks,~,~] = plot_misfit_ks_test_map(dobs,dpred,data_to_compare,u.significance_level,1);
                print_figure(['statistics_compare_',dpred.niter],['ks_test_pval_',num2str(ks.p_all)]); %Save figure
                
            elseif data_main_menu == 11 % cross-plot with a third data response----------------------------------------------
                
                close all
                %This is currently only implemented for loading ModEM data               
                %Load a third data set
                curdir = pwd;               
                [filename,datpath]=uigetfile({'*.dat;*.data', 'ModEM data files (*.dat,*.data)'; '*.*', 'All Files (*.*)'},'Choose a ModEM data file to compare to'); if filename == 0; return; end
                cd(datpath);
                [data_to_compare] = load_data_modem(filename);   
                cd(curdir);
                                
                plot_misfit_cross_plot(dobs,dpred,data_to_compare,1);
                print_figure(['statistics_compare_',dpred.niter],['residuals_cross_plot']); %Save figure
                
            else
                menu3 = 0;


            end


        end

    elseif main_menu == 3 %EDIT MODEL----------------------------------------------------------------

        
        edit_model(model_to_plot,dobs);

              
    elseif main_menu == 4 %EXPORT/EDIT DATA----------------------------------------------------------
        
        while 1
        e_data_menu = menu('','Export Predicted Data as EDI','Export Predicted Data Structure as MAT file','Edit Predicted Data','Quit');
        
        if e_data_menu == 1
            if ~exist('edi_files','dir')
                mkdir('edi_files');
            end
            curdir = pwd;
            cd('edi_files')
            
            for is = 1:dpred.ns
                write_data_edi(dpred,is)        
            end
            
            cd(curdir)
            
        elseif e_data_menu == 2
            
            filename = [basefile,'_pred.mat'];
            save(filename,dpred)
            
        elseif e_data_menu == 3
            
            [dpred] = edit_data(dpred);
            
        else
            break
            
            
        end
        end
            
    elseif main_menu == 5 %RELOAD USER DEFAULTS FILE-----------------------
        clear 'user_defaults.m' % this is needed to load changes saved in editor while S3D is running
        u = user_defaults;
        
        u_path = which('user_defaults.m');       
        disp(['Reloaded user_defaults from ',u_path])
                
    else
        close all
        menu1 = 0;

    end %END MAIN MENU IF

end % END MAIN WHILE LOOP

end %END MAIN

