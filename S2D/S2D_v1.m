function S2D_v1(mode,varargin)
%Function to plot 2D results from Rodi and Mackei (2001) NLCG 2D Inversion.
%
% Apr. 2020: new function detailed_statistics_2D computes misfit for rho
% and phase. the misfit values are usually within 2 decimal places of the
% values in the 2D inversion .rsp file.
%
% 2019: now specify input mode as 'iso' or 'aniso' to try to make this code compatible with iso/aniso inversions
%In order to run this code, you will need two input files in your **current directory**:
%
% 1) The *.par file containing all the inversion settings and station locations
% 2) The *.stn file containing the station names and locations
%
% The code loads the observed data, predicted data, starting model, and
% inversion model into the standard structure formats (dobs, dpred, m, and
% m0)
%
% Potential Bugs or Things to Add:-----------------------------------------
%
% -Need to add static shifts

close all

%------------------User input structures-----------------------------------
%
u = user_defaults;

%--------------------------------------------------------------------------

%
%-----------------------BEGIN LOADING DATA---------------------------------
    
if nargin ==0 || ~strcmp(mode,'iso') && ~strcmp(mode,'aniso')
    error('Must specify input mode as ''iso'' or ''aniso''')
end
    
if isempty(varargin) == 1 %If there is no parameter file specified then get user input

    [startupfile] = uigetfile({'*.par'},'Pick NLCG Inversion Parameter file'); if startupfile == 0; return; end
    [sitefile] = uigetfile({'*.stn'},'Pick NLCG Station file'); if sitefile == 0; return; end

else %If there is a parameter file specified then read the parameter file

    parfile = varargin{1};

    fid = fopen(parfile);
    if fid==-1
        error([parfile,' parameter file not found'])
    else
        reading = textscan(fid,'%s');
        reading = reading{1};
    end
    fclose all;


    startupfile = reading{1}; %Input model file
    sitefile = reading{2}; %Stn file
    
    if length(reading)>2
        site_lat_long_file = reading{3};
    else
        site_lat_long_file = 'none';
    end

end

%------------------LOAD ALL THE DATA-------------------------------------

[par]=load_startup_2DNLCG(mode,startupfile); %load inversion parameter file

if strcmp(mode,'iso')
    m_file = [par.profile_name,'_6_11.mdl'];
    [m0] = load_model_2DNLCG(par.m0_file,sitefile); %load starting model file
    [m] = load_model_2DNLCG(m_file,sitefile); %load inversion model file
    inv_log = load_logfile_2DNLCG([par.profile_name,'_6_11.log']);
else
    m_file = [par.profile_name,'_6_11_xx.mdl'; par.profile_name,'_6_11_yy.mdl'; par.profile_name,'_6_11_zz.mdl'];
    [m0_xx] = load_model_2DNLCG(par.m0_xx_file,sitefile); % load starting model files
    [m0_yy] = load_model_2DNLCG(par.m0_yy_file,sitefile);
    [m0_zz] = load_model_2DNLCG(par.m0_zz_file,sitefile);
    m0 = m0_xx; % combine all 3 models into one m structure
    m0.A_xx = m0.A;
    m0 = rmfield(m0,'A');
    m0.A_yy = m0_yy.A;
    m0.A_zz = m0_zz.A;
    [m_xx] = load_model_2DNLCG(m_file(1,:),sitefile); % load inversion model files
    [m_yy] = load_model_2DNLCG(m_file(2,:),sitefile);
    [m_zz] = load_model_2DNLCG(m_file(3,:),sitefile);    
    m = m_xx;
    m.A_xx = m.A;
    m = rmfield(m,'A');
    m.A_yy = m_yy.A;
    m.A_zz = m_zz.A;
end

dpredfile = [par.profile_name,'_6_11.rsp'];
[dobs, dpred] = load_data_2DNLCG(par,dpredfile,sitefile); %load observed data and predicted data
[dobs, dpred] = get_station_loc_2D(par,dobs,dpred,m);

[s] = detailed_statistics_2D(dobs,dpred);

% dobs.loc = [-dobs.z dobs.y]; % debugging, these are temporary filler
%-----------------DONE LOADING DATA----------------------------------------

%--------------------------------------------------------------------------

%Option to plot either the input model or the inversion model
if u.plot_input_model
    model_to_plot = m0;
    dobs.niter = 'input'; dpred.niter = 'input';
    m0.niter = 'input';
else
    model_to_plot = m; 
    niter = m_file(1,(length(m_file)-6):end-4); 
    
    dobs.niter = niter; dpred.niter = niter;
    m.niter = niter;         
end



%[stats_all,rms_site,rms_freq,rms_comp,rms_sf,rms_fr,residuals] = detailed_statistics(dobs,dpred);

%------------------------BEGINNING MAIN MENU-------------------------------
menu1 = 1;

while menu1 == 1
    
    main_menu = menu('Main Menu','View Model','View Data/Response','Edit Model','Edit/Export Data','Quit');

    if main_menu == 1 %VIEW MODEL------------------------------------------------------------------------------
        menu2 = 1;

        while menu2 == 1
            model_menu = menu('','Vertical Cross-Section','Export Current Cross Section in Lat-Long NOT WORKING','Conductance NOT WORKING','Rho-Depth Curves','Export All Slices in Lat-Long','Return');

            if model_menu==1 %Plot Vertical Cross-Sections
                close all
                set_figure_size(1); hold on
                %Set the ranges to be plotted on the cross-section
                yind = nearestpoint(u.lim_2D(1)*1000,m.cy):nearestpoint(u.lim_2D(2)*1000,m.cy);
                zind = 1:nearestpoint(u.lim_2D(4)*1000,m.cz);
                                              
                if strcmp(mode,'iso')
                    plot_cross_section(m.cy(yind)/1000,m.cz(zind)/1000,m.A(1,yind,zind));  
                    
                    plot(dobs.y/1000,dobs.z/1000,'vk','MarkerFaceColor','k')                
                    title(strrep(['Cross Section. Model ',par.profile_name,' | tau = ',num2str(par.tau),' | rms = ',num2str(s.rms)],'_','\_'));
                    print_figure('vertical_profiles',['Cross_Section_',m.niter]);
                else
                    plot_cross_section(m.cy(yind)/1000,m.cz(zind)/1000,m.A_xx(1,yind,zind));  
                    plot(dobs.y,dobs.z-0.05,'vk','MarkerFaceColor','k')                
                    title(strrep(['Cross Section \rho_{xx} Model ',par.profile_name,' | tau = ',num2str(par.tau),' | aniso tau = ',num2str(par.aniso_tau),' | rms = ',num2str(s.rms)],'_','\_'));
                    print_figure('vertical_profiles',['Cross_Section_rho_xx']);
                    clf
                    
                    plot_cross_section(m.cy(yind)/1000,m.cz(zind)/1000,m.A_yy(1,yind,zind)); 
                    plot(dobs.y,dobs.z-0.05,'vk','MarkerFaceColor','k')                
                    title(strrep(['Cross Section \rho_{yy} Model ',par.profile_name,' | tau = ',num2str(par.tau),' | aniso tau = ',num2str(par.aniso_tau),' | rms = ',num2str(s.rms)],'_','\_'));
                    print_figure('vertical_profiles',['Cross_Section_rho_yy']);
                    clf
                    
                    plot_cross_section(m.cy(yind)/1000,m.cz(zind)/1000,m.A_zz(1,yind,zind));  
                    plot(dobs.y,dobs.z-0.05,'vk','MarkerFaceColor','k')                
                    title(strrep(['Cross Section \rho_{zz} Model ',par.profile_name,' | tau = ',num2str(par.tau),' | aniso tau = ',num2str(par.aniso_tau),' | rms = ',num2str(s.rms)],'_','\_'));
                    print_figure('vertical_profiles',['Cross_Section_rho_zz']);
                    close all
                    
                end
                
            elseif model_menu == 2 %Save currently open cross section
                g=groot;
                if isempty(g.Children)
                    disp('In order to export a N-S, E-W, or diagonal section, you must first choose to plot one')
                else
                    %Write out the section in latitude and
                    %longtiude.
                    write_section(x,y,z,rho,dobs);

                    %Plot the section that you just output (for
                    %debugging purposes)
                    %plot_section_lat_long
                end
            elseif model_menu == 2 %Plot Conductance-----------------------
                
                plot_conductance_2d(model_to_plot,u,dobs)
                print_figure(['conductance',u.num],'conductance_map') 
                    
            elseif model_menu == 3 %Plot resistivity and depth curves------
                
                 rho_z = menu('','View Rho-z Curve Beneath Stations','Pick Model Cell to plot Rho-z Curve','Back');

                 if rho_z ==1
                    plot_rho_z_station(model_to_plot,dobs)
                    
                 elseif rho_z == 2

                    plot_rho_z_model(model_to_plot,dobs)
                    
                 end
                  
            elseif model_menu == 4 %Export slices and cube in lat-long
                
                write_slices_cube_lat_long(model_to_plot);
                
            else
                menu2 = 0;

            end

        end

    elseif main_menu == 2 %VIEW DATA------------------------------------------------------------------------------
        menu3 = 1;

        while menu3 == 1
            % 5 options
            data_main_menu = menu('Data/Response Menu','Plot Data Curves','Plot Data Pseudo','Plot Data Map (NOT WORKING)','Misfit Statistics','Return');
            
            if data_main_menu == 1 % Plot data as curves-------------------
                
                menu5 = menu('Choose Data','Rho / Phase','Impedance','Tipper','Return');
            
                if menu5 == 1 || menu5 == 2 || menu5 == 3 %If you want to plot data as curves                                                    
                    plot_data(dobs,menu5,dpred,s)
                else
                    break
                end
                    
            elseif data_main_menu == 2 % Plot data as pseudo---------------
                
                menu5 = menu('Choose Data','Rho / Phase','Tipper','Return');
                
                if menu5 == 1 % rho / pha               
                    plot_rho_pha_pseudo(dobs)
                    plot_rho_pha_pseudo(dpred)                   
                elseif menu5 == 2 % tipper                  
                    plot_tipper_pseudo(dobs)
                    plot_tipper_pseudo(dpred)                   
                else
                    break                  
                end
                   
            elseif data_main_menu == 3  %Plot data in Map View-------------------------------
                    
                menu5 = menu('Choose Data','Rho / Phase','Tipper','Return');
                
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

                        if menu5 == 1 %Plot interpolated rho and phase in map view
                            plot_rho_pha_map(dobs,ip)
                            plot_rho_pha_map(dpred,ip)
%                             plot_percent_difference_map(dobs,dpred,ip,u)
                        elseif menu5 == 2 %Plot induction vectors in map view
                            close all
                            plot_induction_vector_map(dobs,ip,'r')
                            plot_induction_vector_map(dpred,ip,'k')
                            figure_name = 'induction_vector';
                            manual_legend('Observed Data','-r','Modelled Data','-k');
                            title(['Induction Vectors for T = ',num2str(dobs.T(ip))])
                            print_figure(['compare_responses_map_',u.num],[figure_name,'_',num2str(dobs.T(ip))]); %Save figure
                        else
                            break
                        end
                    end

            elseif data_main_menu == 4 %Misfit Statistics----------------------------------------------

                misfit_menu = menu('','Residual Histograms','Misfit By Period','Misfit By Inversion Iteration');


                if misfit_menu == 1 %Residual histograms-------------------

                    plot_misfit_residual_histograms(s)
                    print_figure(['misfit_stats_',u.num],'Residuals'); %Save figure

                elseif misfit_menu == 2 %Misfit by period------------------

                    plot_misfit_rms_by_period(dobs,dpred,s)
                    print_figure(['misfit_stats_',u.num],'rms_by_period_by_component'); %Save figure
                    
                elseif misfit_menu == 3 %Misfit convergence----------------
                    plot_misfit_convergence(inv_log)
                    print_figure(['misfit_stats_',u.num],'inversion_convergence'); %Save figure
                end


            else
                menu3 = 0;


            end


        end

    elseif main_menu == 3 %EDIT MODEL----------------------------------------------------------------

        %edit_model(model_to_plot)

            
        
    elseif main_menu == 4 %EXPORT/EDIT DATA
        
        
        e_data_menu = menu('','Export Predicted Data as EDI','Export Predicted Data as MAT file','Edit Predicted Data','Quit');
        
        if e_data_menu == 1
            
            %write_edi(dpred)
            
        elseif e_data_menu == 2
            
            %write_mat(dpred)
            
        elseif e_data_menu == 3
            
            %edit_data(dpred)
            
            
        end
            


    else

        menu1 = 0;

    end %END MAIN MENU IF

end % END MAIN WHILE LOOP

close all
end %END MAIN

