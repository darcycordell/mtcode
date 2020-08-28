function M3D_v1
% Re-write M3 without GUIDE
%
% This function is the main GUI program used to build MT model meshes and
% output data and models in formats that can be read by external MT inversion
% scripts such as ModEM and WSINV3D
%
% The primary input needed for this program is a '*.mat' file which
% contains all the EDI data for a given MT survey. This mat file can be created
% user interpolate_edi.m. See interpolate_edi.m for more information.
% New Oct. 2019 - You can import a ModEM format data file instead of a .mat file.
%
% Once the data file is loaded, the user can specify mesh dimensions,
% add topography, view data, output data etc.
%
%
% (C) 2020 Unsworth Research Group (University of Alberta, Edmonton, Canada)

M3D_v1_OpeningFcn
function M3D_v1_OpeningFcn
    
    clear global
    clc    
    
    H = [];
%     f = figure('name','M3D_v1','numbertitle','off','units','pixels','Position',[200 200 1122 768]);
    f = figure('name','M3D_v1','numbertitle','off','units','normalized','Position',[0.1 0.1 0.9 0.9],'outerposition',[0.1 0.1 0.9 0.9]);
    mp = 1./[1122 768 1122 768]; % this is the original size of the GUI. It is used as a factor to normalize all components to the figure size
    delete undo*
    if isfield(H,'undo')
        H = rmfield(H,'undo');
        H.undo=0;
    else
        H.undo=0;
    end
%     set(f,'Toolbar','figure');  % Display the standard toolbar
    disp('THIS PROGRAM PREPARES 3D MESH GRID AND DATA FILES FOR MODEM AND WSINV3DMT')
    disp('Based on M3DET by Ersan Turkoglu (2006)')
    disp('University of Alberta, 2018')
   
    axes1 = axes('parent',f,'units','normalized','Position',[42 68 644 586].*mp);
    [x,y,z] = meshgrid(-2:.2:2, -2:.25:2, -2:.16:2);
    v = x .* exp(-x.^2 - y.^2 - z.^2);
    slice(x,y,z,v,[-1.2 .8 2],2,[-2 -.2])
    axis equal; axis tight; clear x y z
    xlb='Easting (km)'; ylb='Northing (km)'; zlb='Depth (km)';
    xlabel(xlb);ylabel(ylb);zlabel(zlb)
    axis off; % Turn off all axis labeling   
    set(axes1,'HandleVisibility','OFF');
%%  LOAD PANEL
% This panel involves loading data, loading models, and the M3D parameter file.

    panel_p = [734 644 380 110];
    comp_p = [panel_p(3) panel_p(4) panel_p(3) panel_p(4)];
    uipanel_load_menu = uipanel('Title','1. LOAD FILES','FontSize',8,'units','normalized','Position',panel_p.*mp);
        
        uicontrol('Parent',uipanel_load_menu,'String','Load .mat Data File','units','normalized','Position',[10 72 170 22]./comp_p,'Callback',@M3_load_mat_new);    
        uicontrol('Parent',uipanel_load_menu,'String','Load .mat Data File (Legacy)','units','normalized','Position',[10 49 170 22]./comp_p,'Callback',@M3_load_mat_legacy);
        uicontrol('Parent',uipanel_load_menu,'String','Load ModEM Data File','units','normalized','Position',[10 26 170 22]./comp_p,'Callback',@M3_load_modem_data);
        
        uicontrol('Parent',uipanel_load_menu,'String','Load wsinv3dmt Model','units','normalized','Position',[200 72 170 22]./comp_p,'Callback',@M3_load_ws_mod);
        uicontrol('Parent',uipanel_load_menu,'String','Load Winglink Model','units','normalized','Position',[200 49 170 22]./comp_p,'Callback',@M3_load_winglink_mod);
        uicontrol('Parent',uipanel_load_menu,'String','Load ModEM Model','units','normalized','Position',[200 26 170 22]./comp_p,'Callback',@M3_load_modem_mod);
        
        uicontrol('Parent',uipanel_load_menu,'String','Load M3D Parameter File','units','normalized','Position',[128 3 128 22]./comp_p,'Callback',@M3_load_param);
        
%%  MAKE MODEL PANEL
% This panel involves generating a model.

    panel_p = [734 472 380 170];
    comp_p = [panel_p(3) panel_p(4) panel_p(3) panel_p(4)];
    uipanel_model = uipanel('Title','2. MAKE MODEL','FontSize',8,'units','normalized','Position',panel_p.*mp);
       
        uicontrol('Parent',uipanel_model,'style','text','string','Cell width (km)','units','normalized','Position',[8 132 80 16]./comp_p);
        uicontrol('Parent',uipanel_model,'style','text','string','Number of padding cells','units','normalized','Position',[8 108 80 16]./comp_p);
        uicontrol('Parent',uipanel_model,'style','text','string','Padding factor','units','normalized','Position',[8 84 80 16]./comp_p);
        uicontrol('Parent',uipanel_model,'style','text','string','X','units','normalized','Position',[111 155 16 16]./comp_p);
        uicontrol('Parent',uipanel_model,'style','text','string','Y','units','normalized','Position',[141 155 16 16]./comp_p);
        
%         uicontrol('Parent',uipanel_model,'style','text','string','Spacing in Y (km)','units','normalized','Position',[188 132 116 16]./comp_p);
%         uicontrol('Parent',uipanel_model,'style','text','string','Number of padding in Y','units','normalized','Position',[188 108 116 16]./comp_p);
%         uicontrol('Parent',uipanel_model,'style','text','string','Increase by in Y','units','normalized','Position',[188 84 116 16]./comp_p);
        
        uicontrol('Parent',uipanel_model,'style','text','string','First thickness (m)','units','normalized','Position',[188 132 116 16]./comp_p);
        uicontrol('Parent',uipanel_model,'style','text','string','Increase thickness by','units','normalized','Position',[188 108 116 16]./comp_p);
        uicontrol('Parent',uipanel_model,'style','text','string','Number of Layers','units','normalized','Position',[188 84 116 16]./comp_p);            
        
        H.xspacing = uicontrol('Parent',uipanel_model,'style','edit','string','5','backgroundcolor','white','tooltip',sprintf('Enter cell size in core mesh.'),'units','normalized','Position',[107 132 24 22]./comp_p,'Callback',@xspacing_Callback);
        H.yspacing = uicontrol('Parent',uipanel_model,'style','edit','string','5','backgroundcolor','white','units','normalized','Position',[137 132 24 22]./comp_p,'Callback',@yspacing_Callback); 
        H.moutcX = uicontrol('Parent',uipanel_model,'style','edit','string','8','backgroundcolor','white','units','normalized','Position',[107 108 24 22]./comp_p,'Callback',@moutcX_Callback);
        H.moutcY = uicontrol('Parent',uipanel_model,'style','edit','string','8','backgroundcolor','white','units','normalized','Position',[137 108 24 22]./comp_p,'Callback',@moutcY_Callback);
        H.Xinc = uicontrol('Parent',uipanel_model,'style','edit','string','1.3','backgroundcolor','white','units','normalized','Position',[107 84 24 22]./comp_p,'Callback',@Xinc_Callback);                         
        H.Yinc = uicontrol('Parent',uipanel_model,'style','edit','string','1.3','backgroundcolor','white','units','normalized','Position',[137 84 24 22]./comp_p,'Callback',@Yinc_Callback);
        
        H.firstthick = uicontrol('Parent',uipanel_model,'style','edit','string','100','backgroundcolor','white','units','normalized','Position',[307 132 30 22]./comp_p,'Callback',@firstthick_Callback);
        H.Zinc = uicontrol('Parent',uipanel_model,'style','edit','string','1.3','backgroundcolor','white','units','normalized','Position',[307 108 30 22]./comp_p,'Callback',@Zinc_Callback);       
        H.nl = uicontrol('Parent',uipanel_model,'style','edit','string','33','backgroundcolor','white','units','normalized','Position',[307 84 30 22]./comp_p,'Callback',@nl_Callback);       
                         
        uicontrol('Parent',uipanel_model,'style','text','string','Model Resistivity (Ohm-m)','units','normalized','Position',[5 50 100 21]./comp_p)
        H.res = uicontrol('Parent',uipanel_model,'style','edit','string','100','units','normalized','Position',[100 55 31 21]./comp_p,'Callback',@res_Callback);  
        
        uicontrol('Parent',uipanel_model,'style','pushbutton','string','Generate Model','units','normalized','Position',[20 26 88 22]./comp_p,'Callback',@M3_generate_model);
        uicontrol('Parent',uipanel_model,'style','pushbutton','string','View Model','units','normalized','Position',[20 2 88 22]./comp_p,'Callback',@M3_view_model);

        subpanel_rot = uipanel('parent',uipanel_model,'Title','Rotation','FontSize',8,'units','normalized','Position',[0.4 0 0.3 0.45]);
            uicontrol('Parent',subpanel_rot,'style','text','string','Rotation Angle','tooltip','Enter strike angle, positive = clockwise from north','units','normalized','Position',[0 0.575 0.8 0.25])
            H.rotang = uicontrol('Parent',subpanel_rot,'style','edit','string','0','tooltip','Enter strike angle, positive = clockwise from north','units','normalized','Position',[0.75 0.55 0.2 0.3],'Callback',@rotang_Callback);  
            uicontrol('Parent',subpanel_rot,'style','pushbutton','string','ROTATE STATIONS & DATA','tooltip',sprintf('Enter angle, positive = clockwise from north\nData are rotated to this angle\nMesh x-direction (North) will be aligned to this angle'),'units','normalized','Position',[0.04 0.1 0.92 0.35],'Callback',@M3_set_rot);

        
        subpanel_topo = uipanel('parent',uipanel_model,'Title','Add Topography','FontSize',8,'units','normalized','Position',[0.7 0 0.3 0.45]);
            uicontrol('Parent',subpanel_topo,'style','text','string','Air cell thickness (m)','units','normalized','Position',[0 0.6 0.8 0.2]);
            H.airspacing = uicontrol('Parent',subpanel_topo,'style','edit','string','20','backgroundcolor','white','units','normalized','Position',[0.75 0.55 0.2 0.3],'Callback',@airspacing_Callback);
            uicontrol('Parent',subpanel_topo,'style','pushbutton','string','Add Topo','units','normalized','Position',[0.1 0.1 0.8 0.35],'Callback',@M3_add_topo);
        
        %%  EDIT MODEL PANEL     
% This panel involves editing the model.  

    panel_p = [734 355 380 116];
    comp_p = [panel_p(3) panel_p(4) panel_p(3) panel_p(4)];
    uipanel_edit_model = uipanel('Title','3. EDIT MODEL','FontSize',8,'units','normalized','Position',panel_p.*mp);
       
        uicontrol('Parent',uipanel_edit_model,'style','pushbutton','string','Add Row','units','normalized','Position',[6 74 78 22]./comp_p,'Callback',@M3_add_row)
        uicontrol('Parent',uipanel_edit_model,'style','pushbutton','string','Add Column','units','normalized','Position',[90 74 78 22]./comp_p,'Callback',@M3_add_col)
        uicontrol('Parent',uipanel_edit_model,'style','pushbutton','string','Delete Row','units','normalized','Position',[6 50 78 22]./comp_p,'Callback',@M3_del_row)
        uicontrol('Parent',uipanel_edit_model,'style','pushbutton','string','Delete Column','units','normalized','Position',[90 50 78 22]./comp_p,'Callback',@M3_del_col)

        uicontrol('Parent',uipanel_edit_model,'style','pushbutton','string','Prev. Layer','units','normalized','Position',[200 74 78 22]./comp_p,'Callback',@M3_prev_layer)
        uicontrol('Parent',uipanel_edit_model,'style','pushbutton','string','Next Layer','units','normalized','Position',[284 74 78 22]./comp_p,'Callback',@M3_next_layer)
        uicontrol('Parent',uipanel_edit_model,'style','pushbutton','string','Jump to Layer','units','normalized','Position',[200 50 78 22]./comp_p,'Callback',@M3_jump_to_layer)
        H.lay_num = uicontrol('Parent',uipanel_edit_model,'style','edit','string','1','backgroundcolor','white','units','normalized','Position',[284 50 30 22]./comp_p,'Callback',@lay_num_Callback);
        
        uicontrol('Parent',uipanel_edit_model,'style','pushbutton','string','Smooth Last Air Layer','tooltip',(''),'units','normalized','Position',[6 2 164 22]./comp_p,'Callback',@M3_smooth_air);
        uicontrol('Parent',uipanel_edit_model,'style','pushbutton','string','Delete Layers','tooltip',(''),'units','normalized','Position',[6 26 164 22]./comp_p,'Callback',@M3_delete_layers);
        
%%  OUTPUT DATA PANEL      
% This panel involves selecting data periods and setting error floors.

    panel_p = [734 246 380 108];
    comp_p = [panel_p(3) panel_p(4) panel_p(3) panel_p(4)];
    uipanel_outputdata = uipanel('Title','4. OUTPUT DATA','FontSize',8,'units','normalized','Position',panel_p.*mp);   
    
        H.min_per = uicontrol('Parent',uipanel_outputdata,'style','edit','string','0.001','backgroundcolor','white','tooltip',sprintf('Enter minimum period to be exported.\nNo periods below this value will be exported.'),'units','normalized','Position',[6 54 36 22]./comp_p,'Callback',@min_per_Callback);
        H.max_per = uicontrol('Parent',uipanel_outputdata,'style','edit','string','1000','backgroundcolor','white','tooltip',sprintf('Enter maximum period to be exported.\nNo periods above this value will be exported.'),'units','normalized','Position',[66 54 36 22]./comp_p,'Callback',@max_per_Callback);
        H.per_skip = uicontrol('Parent',uipanel_outputdata,'style','edit','string','1','backgroundcolor','white','tooltip',sprintf('Export every nth period.\nExample: "1" takes all periods, "2" takes every other...'),'units','normalized','Position',[46 54 18 22]./comp_p,'Callback',@per_skip_Callback);
        H.errflr_Zdiag = uicontrol('Parent',uipanel_outputdata,'style','edit','string','10','backgroundcolor','white','units','normalized','Position',[112 54 18 22]./comp_p,'Callback',@errflr_Zdiag_Callback);
        H.errflr_Zodiag = uicontrol('Parent',uipanel_outputdata,'style','edit','string','10','backgroundcolor','white','units','normalized','Position',[112 30 18 22]./comp_p,'Callback',@errflr_Zodiag_Callback);
        H.errflr_type = uicontrol('Parent',uipanel_outputdata,'style','edit','string','1','backgroundcolor','white','tooltip',sprintf('Type 1: All Z components error = %% of sqrt(|Zxy*Zyx|)\nType 2: Zxx error = %% of |Zxy|, and so on'),'units','normalized','Position',[112 6 18 22]./comp_p,'Callback',@errflr_type_Callback);
        H.errflr_tip = uicontrol('Parent',uipanel_outputdata,'style','edit','string','0.02','backgroundcolor','white','units','normalized','Position',[200 54 28 22]./comp_p,'Callback',@errflr_tip_Callback);
        
        uicontrol('Parent',uipanel_outputdata,'style','text','string','Periods','units','normalized','Position',[30 80 50 15]./comp_p);
        uicontrol('Parent',uipanel_outputdata,'style','text','string','% Error Floor','units','normalized','Position',[110 76 70 15]./comp_p);
        uicontrol('Parent',uipanel_outputdata,'style','text','string','Zxx Zyy','units','normalized','Position',[132 49 43 23]./comp_p);
        uicontrol('Parent',uipanel_outputdata,'style','text','string','Zxy Zyx','units','normalized','Position',[132 25 43 23]./comp_p);
        uicontrol('Parent',uipanel_outputdata,'style','text','string','Err flr type','tooltip',sprintf('Type 1: All Z components error = %% of sqrt(|Zxy*Zyx|)\nType 2: Zxx error = %% of |Zxy|, and so on'),'units','normalized','Position',[130 1 43 23]./comp_p);
        uicontrol('Parent',uipanel_outputdata,'style','text','string','Tipper err','units','normalized','Position',[230 49 54 23]./comp_p);
        
        H.resp = uicontrol('Parent',uipanel_outputdata,'style','popupmenu','string',{'Full Tensor (8)';'FT + Tip (12)';'Off-Diag Z Only (4)';'Tipper Only (4)'},'units','normalized','Position',[6 20 96 20]./comp_p,'Callback',@num_resp_Callback); 
        H.conj_Z = uicontrol('Parent',uipanel_outputdata,'style','checkbox','string','Conj. Z','tooltip',sprintf('WARNING: matfile time dep. is e^+iwt,\n wsinv3dmt needs e^-iwt, and ModEM can use either one.'),'units','normalized','Position',[200 28 60 23]./comp_p,'Callback',@conj_Z_Callback);
        H.conj_tip = uicontrol('Parent',uipanel_outputdata,'style','checkbox','string','Conj. Tip','units','normalized','Position',[200 8 60 23]./comp_p,'Callback',@conj_tip_Callback);
        

        uicontrol('Parent',uipanel_outputdata,'style','pushbutton','string','View Data','units','normalized','Position',[290 12 76 76]./comp_p,'Callback',@M3_view_data);
%%  SAVE PANEL
% This panel involves saving data, model, and convariance files.

    panel_p = [834 100 280 144];
    comp_p = [panel_p(3) panel_p(4) panel_p(3) panel_p(4)];
    uipanel_save_menu = uipanel('Title','5. SAVE FILES','FontSize',8,'units','normalized','Position',panel_p.*mp);
    
    uipanel_save_wsinv = uipanel('parent',uipanel_save_menu,'title','WSINV3DMT','FontSize',8,'units','normalized','Position',[0.52 0.2 0.46 0.8]);
            
        uicontrol('Parent',uipanel_save_wsinv,'style','pushbutton','string','Save wsinv3dmt Data','units','normalized','Position',[0.1 0.7 0.8 0.25],'Callback',@M3_save_ws_data);      
        uicontrol('Parent',uipanel_save_wsinv,'style','pushbutton','string','Save wsinv3dmt Model','units','normalized','Position',[0.1 0.4 0.8 0.25],'Callback',@M3_save_ws_mod);
        uicontrol('Parent',uipanel_save_wsinv,'String','Save WS Startup File','units','normalized','Position',[0.1 0.1 0.8 0.25],'Callback',@M3_make_ws_startup);
     
    uipanel_save_modem = uipanel('parent',uipanel_save_menu,'title','ModEM','FontSize',8,'units','normalized','Position',[0.02 0.2 0.46 0.8]);  
    
        uicontrol('Parent',uipanel_save_modem,'style','pushbutton','string','Save ModEM Data','units','normalized','Position',[0.1 0.7 0.8 0.25],'Callback',@M3_save_modem_data);
        uicontrol('Parent',uipanel_save_modem,'style','pushbutton','string','Save ModEM Model','units','normalized','Position',[0.1 0.4 0.8 0.25],'Callback',@M3_save_modem_mod);
        uicontrol('Parent',uipanel_save_modem,'String','Save ModEM cov File','units','normalized','Position',[0.1 0.1 0.8 0.25],'Callback',@M3_make_modem_cov);
        
        uicontrol('Parent',uipanel_save_menu,'String','Save M3D Parameter File','units','normalized','Position',[72 4 140 22]./comp_p,'Callback',@M3_make_param);                  
        
%%  OCEAN PANEL
% This panel is not currently used and may be removed in the future.

    uipanel_ocean_model = uipanel('Title','OCEAN MENU','FontSize',8,'units','normalized','Position',[596 667 128 89].*mp);
    
    H.ocean_res = uicontrol('Parent',uipanel_ocean_model,'style','edit','string','0.3','backgroundcolor','white','Position',[10 26 22 20],'Callback',@ocean_res_Callback);
    H.fix_ocean = uicontrol('Parent',uipanel_ocean_model,'style','checkbox','string','Fix Ocean','Position',[10 6 68 17],'Callback',@fix_ocean_Callback);
        
    uicontrol('Parent',uipanel_ocean_model,'style','text','string','Ocean Resistivity','Position',[35 28 85 16]);

%% VIEW OPTIONS PANEL

    panel_p = [834 20 280 78];
    comp_p = [panel_p(3) panel_p(4) panel_p(3) panel_p(4)];
    uipanel_view_options = uipanel('Title','VIEW OPTIONS','FontSize',8,'units','normalized','Position',panel_p.*mp);
        
        uicontrol('Parent',uipanel_view_options,'style','text','string','Resistivity','units','normalized','Position',[0.1 0.7 0.2 0.2])
        uicontrol('Parent',uipanel_view_options,'style','text','string','Min','units','normalized','Position',[0.08 0.45 0.1 0.2])
        uicontrol('Parent',uipanel_view_options,'style','text','string','Max','units','normalized','Position',[0.22 0.45 0.1 0.2])   
        H.min_res = uicontrol('Parent',uipanel_view_options,'style','edit','string','1','backgroundcolor','white','units','normalized','Position',[0.08 0.1 0.1 0.3],'Callback',@M3_min_res);
        H.max_res = uicontrol('Parent',uipanel_view_options,'style','edit','string','1000','backgroundcolor','white','units','normalized','Position',[0.22 0.1 0.1 0.3],'Callback',@M3_max_res);        
          

%         uicontrol('Parent',uipanel_view_options,'style','text','string','Plot Buffer (km)','units','normalized','Position',[0.5 0.65 0.2 0.2])
%         stat_tol = uicontrol('Parent',uipanel_view_options,'style','edit','string','4','units','normalized','Position',[0.75 0.65 0.1 0.3],'Callback',@stat_tol_Callback);     
           
        H.gridlines = uicontrol('Parent',uipanel_view_options,'style','checkbox','string','No Grid','units','normalized','Position',[0.5 0.45 0.2 0.2],'Callback',@gridlines_Callback);
        H.station_names = uicontrol('Parent',uipanel_view_options,'style','checkbox','string','Station Names','units','normalized','Position',[0.5 0.1 0.4 0.2],'Callback',@station_names_Callback);
    
    %%  MISC UICONTROL COMPONENTS
% Note: the positions of these components is relative to the figure size

    H.model_text = uicontrol('style','text','string','Mesh info','FontSize',10,'units','normalized','Position',[120 685 280 80].*mp);
    H.data_text = uicontrol('style','text','string','Number of periods = ','FontSize',10,'units','normalized','Position',[405 700 180 60].*mp);
       
    uicontrol('String','Restart','units','normalized','Position',[12 740 60 22].*mp,'Callback',@restart_button_Callback);
    uicontrol('String','Exit','units','normalized','Position',[12 716 60 22].*mp,'Callback',@exit_Callback);

%%   
    % Update H structure, provide default values
    H.axes1 = axes1;
%     H.stat_tol = get(stat_tol,'string');
    
    H.NR = 8;
    H.per = [];
      
    H.numair = 0; % default zero
      
    guidata(f, H)
end % end opening function

%% Callback Functions Below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function stat_tol_Callback(hObject,~) % IN PROGRESS does nothing
    H=guidata(gcbo);
    H.stat_tol = get(hObject,'string');    
    guidata(hObject, H);    
end

function xspacing_Callback(hObject,~)
    H=guidata(gcbo); 
    set(H.xspacing,'String',get(hObject,'String'))
    guidata(hObject, H); 
end

function yspacing_Callback(hObject,~)
    H=guidata(gcbo);
    set(H.yspacing,'String',get(hObject,'String')) 
    guidata(hObject, H); 
end

function moutcX_Callback(hObject,~)
    H=guidata(gcbo);
    set(H.moutcX,'String',get(hObject,'String')) 
    guidata(hObject, H); 
end

function moutcY_Callback(hObject,~)
    H=guidata(gcbo);
    set(H.moutcY,'String',get(hObject,'String'))  
    guidata(hObject, H); 
end

function firstthick_Callback(hObject,~)
    H=guidata(gcbo);
    set(H.firstthick,'String',get(hObject,'String'))  
    guidata(hObject, H); 
end

function nl_Callback(hObject,~)
    H=guidata(gcbo);
    set(H.nl,'String',get(hObject,'String'))
    guidata(hObject, H); 
end

function Xinc_Callback(hObject,~)
    H=guidata(gcbo);
    set(H.Xinc,'String',get(hObject,'String'))  
    guidata(hObject, H); 
end

function Yinc_Callback(hObject,~)
    H=guidata(gcbo);
    set(H.Yinc,'String',get(hObject,'String')) 
    guidata(hObject, H); 
end

function Zinc_Callback(hObject,~)
    H=guidata(gcbo);
    set(H.Zinc,'String',get(hObject,'String'))   
    guidata(hObject, H); 
end

function airspacing_Callback(hObject,~)
    H=guidata(gcbo);
    set(H.airspacing,'String',get(hObject,'String'))   
    guidata(hObject, H); 
end

function lay_num_Callback(hObject,~)
    H=guidata(gcbo); 
    if str2double(get(hObject,'String')) < 1 || str2double(get(hObject,'String')) > (H.nz - 1)
        warndlg('Enter a valid layer number.')
        set(H.lay_num,'String','1')
    else
        set(H.lay_num,'String',get(hObject,'String'))
    end
    guidata(hObject, H); 
end

%%

function min_per_Callback(hObject,~)    
    H=guidata(gcbo);
    set(H.min_per,'String',get(hObject,'String'))
    guidata(hObject, H);
end

function max_per_Callback(hObject,~)    
    H=guidata(gcbo);
    set(H.max_per,'String',get(hObject,'String'))
    guidata(hObject, H);
end

function per_skip_Callback(hObject,~)    
    H=guidata(gcbo);
    if str2double(get(hObject,'String')) == 0
        warndlg('Warning, enter a value greater than zero.')
        set(H.per_skip,'string','1') 
    else
        set(H.per_skip,'String',get(hObject,'String'))
    end
    guidata(hObject, H);
end

function errflr_Zdiag_Callback(hObject,~)
    H=guidata(gcbo);
    set(H.errflr_Zdiag,'String',get(hObject,'String'))
    guidata(hObject, H);
end

function errflr_Zodiag_Callback(hObject,~)
    H=guidata(gcbo);
    set(H.errflr_Zodiag,'String',get(hObject,'String'))
    guidata(hObject, H);
end

function errflr_type_Callback(hObject,~)
    H=guidata(gcbo);
    set(H.errflr_type,'String',get(hObject,'String'))
    guidata(hObject, H);
end

function errflr_tip_Callback(hObject,~)
    H=guidata(gcbo);
    set(H.errflr_tip,'String',get(hObject,'String'))
    guidata(hObject, H);
end

function num_resp_Callback(hObject,~)
    H=guidata(gcbo);
    value = get(hObject,'value'); % check which popup selection is made
    % Note two variables: resp is a number for each separate option, NR is the number of responses
    if value == 1 % default, full impedance tensor (nr = 8)
        set(H.resp,'Value',get(hObject,'Value')); H.NR = 8;
    elseif value == 2 % full tensor and tip (nr = 12)
        set(H.resp,'Value',get(hObject,'Value')); H.NR = 12;
    elseif value == 3 % off-diag Z only (nr = 4)
        set(H.resp,'Value',get(hObject,'Value')); H.NR = 4;
    elseif value == 4 % tipper only (nr = 4)
        set(H.resp,'Value',get(hObject,'Value')); H.NR = 4;
    end
    guidata(hObject, H);
end

function conj_Z_Callback(hObject,~)
    H=guidata(gcbo);   
    set(H.conj_Z,'Value',get(hObject,'Value')) % see if box is unchecked (0) or checked (1)
    guidata(hObject, H);
end

function conj_tip_Callback(hObject,~)
    H=guidata(gcbo); 
    set(H.conj_tip,'Value',get(hObject,'Value')) % see if box is unchecked (0) or checked (1)
    guidata(hObject, H);
end

function ocean_res_Callback(hObject,~)
    H=guidata(gcbo);
    set(H.ocean_res,'String',get(hObject,'String')) 
    guidata(hObject, H); 
end

function fix_ocean_Callback(hObject,~)
    H=guidata(gcbo);
    set(H.fix_ocean,'Value',get(hObject,'Value')) 
    guidata(hObject, H); 
end

function res_Callback(hObject,~) % IN PROGRESS, currently does nothing
    H=guidata(gcbo);
    set(H.res,'String',get(hObject,'String'))
    guidata(hObject, H); 
end

function rotang_Callback(hObject,~)
    H=guidata(gcbo);
    set(H.rotang,'String',get(hObject,'String'))
    guidata(hObject, H); 
end

%%

function gridlines_Callback(hObject,~)
    H=guidata(gcbo);
    if isfield(H,'m')
        if get(hObject,'Value') % see if box is unchecked (0) or checked (1)
            shading flat
        else
            shading faceted
        end
    else
        warndlg('Create or import a model first.')
        set(hObject,'value',0)
    end
    guidata(hObject, H);
end

function station_names_Callback(hObject,~)
    H=guidata(gcbo);
    if isfield(H,'d')
        if get(hObject,'Value')
            H.tx = text(H.d.y./1000,H.d.x./1000,H.d.site,'fontsize',11,'FontWeight','bold');
%             H.axes2 = axes('Visible','on','hittest','off'); % Invisible axes
%             linkprop([H.axes1 H.axes2],{'CameraPosition' 'XLim' 'YLim' 'ZLim' 'Position'}); % The axes should stay aligned
%             set(H.tx,'Parent',H.axes2); % Put the text in the invisible Axes
%             H.tx = text(H.axes1,H.x,H.y,H.site,'fontsize',11,'FontWeight','bold');
            
            
%             set(H.tx,'handlevisibility','callback')
            set(H.tx,'Clipping','on');
        else
            set(H.tx,'visible','off')
        end
    else
        warndlg('Import data first.')
        set(hObject,'value',0)
    end
    guidata(hObject, H);
end

function restart_button_Callback(~,~)
    clear all; close all; M3D_v1
    delete undo*
end

function exit_Callback(~,~)
    close all; clear all
    delete undo*
end

end % end M3