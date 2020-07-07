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
% Once the mat file is loaded, then the user can specifiy mesh dimensions,
% add topography, view data, output data etc.
%
M3D_v1_OpeningFcn
function M3D_v1_OpeningFcn
    
    clear global
    clc    
    
    H = [];
    f = figure('name','M3D_v1','numbertitle','off','units','pixels','Position',[200 200 1122 768]);
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
   
    axes1 = axes('parent',f,'units','pixels','Position',[50 68 674 586]);
    [x,y,z] = meshgrid(-2:.2:2, -2:.25:2, -2:.16:2);
    v = x .* exp(-x.^2 - y.^2 - z.^2);
    slice(x,y,z,v,[-1.2 .8 2],2,[-2 -.2])
    axis equal; axis tight; clear x y z
    xlb='Easting (km)'; ylb='Northing (km)'; zlb='Depth (km)';
    xlabel(xlb);ylabel(ylb);zlabel(zlb)
    axis off; % Turn off all axis labeling   
    set(axes1,'HandleVisibility','OFF');

    uipanel_menu = uipanel('Title','Menu','FontSize',8,'units','pixels','Position',[734 426 190 328]);
    % This panel involves loading data, loading models, and exporting misc. inversion files.
    
        uicontrol('Parent',uipanel_menu,'String','Load Data from .mat file','Position',[10 290 170 22],'Callback',@M3_load_mat_new);    
        uicontrol('Parent',uipanel_menu,'String','Load Data from .mat file (Legacy)','Position',[10 267 170 22],'Callback',@M3_load_mat_legacy);
        uicontrol('Parent',uipanel_menu,'String','Load ModEM Data File','Position',[10 244 170 22],'Callback',@M3_load_modem_data);
        uicontrol('Parent',uipanel_menu,'String','Load wsinv3dmt Model','Position',[17 167 140 22],'Callback',@M3_load_ws_mod);
        uicontrol('Parent',uipanel_menu,'String','Load Winglink Model','Position',[17 144 140 22],'Callback',@M3_load_winglink_mod);
        uicontrol('Parent',uipanel_menu,'String','Load ModEM model','Position',[17 121 140 22],'Callback',@M3_load_modem_mod);
        uicontrol('Parent',uipanel_menu,'String','Make ModEM cov File','Position',[27 98 120 22],'Callback',@M3_make_modem_cov);
        uicontrol('Parent',uipanel_menu,'String','Make WS Startup File','Position',[27 75 120 22],'Callback',@M3_make_ws_startup);
        uicontrol('Parent',uipanel_menu,'String','Load M3 Parameter File','Position',[27 52 120 22],'Callback',@M3_load_param);
        uicontrol('Parent',uipanel_menu,'String','Make M3 Parameter File','Position',[27 29 120 22],'Callback',@M3_make_param);
        uicontrol('Parent',uipanel_menu,'String','Restart','Position',[27 6 60 22],'Callback',@restart_button_Callback);
        uicontrol('Parent',uipanel_menu,'String','Exit','Position',[87 6 60 22],'Callback',@exit_Callback);
       
    uipanel_model = uipanel('Title','Model','FontSize',8,'units','pixels','Position',[928 401 178 353]);
    % This panel involves generating and saving models.
    
        uicontrol('Parent',uipanel_model,'style','text','string','Spacing in X (km)','Position',[8 306 88 16]);
        uicontrol('Parent',uipanel_model,'style','text','string','Spacing in Y (km)','Position',[8 282 89 16]);
        uicontrol('Parent',uipanel_model,'style','text','string','Mesh out the core in X','Position',[8 258 113 16]);
        uicontrol('Parent',uipanel_model,'style','text','string','Increase by in X','Position',[8 234 81 16]);
        uicontrol('Parent',uipanel_model,'style','text','string','Mesh out the core in Y','Position',[8 210 113 16]);
        uicontrol('Parent',uipanel_model,'style','text','string','Increase by in Y','Position',[8 186 86 16]);
        uicontrol('Parent',uipanel_model,'style','text','string','First thickness (m)','Position',[8 162 93 16]);
        uicontrol('Parent',uipanel_model,'style','text','string','Number of Layers','Position',[8 138 89 16]);
        uicontrol('Parent',uipanel_model,'style','text','string','Increase thickness by','Position',[8 114 111 16]);
        uicontrol('Parent',uipanel_model,'style','text','string','Air cell thickness (m)','Position',[8 90 111 16]);
        
        H.xspacing = uicontrol('Parent',uipanel_model,'style','edit','string','5','backgroundcolor','white','tooltip',sprintf('Enter cell size in core mesh.'),'Position',[127 303 30 22],'Callback',@xspacing_Callback);
        H.yspacing = uicontrol('Parent',uipanel_model,'style','edit','string','5','backgroundcolor','white','Position',[127 279 30 22],'Callback',@yspacing_Callback);
        H.moutcX = uicontrol('Parent',uipanel_model,'style','edit','string','8','backgroundcolor','white','Position',[127 255 30 22],'Callback',@moutcX_Callback);
        H.Xinc = uicontrol('Parent',uipanel_model,'style','edit','string','1.3','backgroundcolor','white','Position',[127 231 30 22],'Callback',@Xinc_Callback);
        H.moutcY = uicontrol('Parent',uipanel_model,'style','edit','string','8','backgroundcolor','white','Position',[127 207 30 22],'Callback',@moutcY_Callback);
        H.Yinc = uicontrol('Parent',uipanel_model,'style','edit','string','1.3','backgroundcolor','white','Position',[127 183 30 22],'Callback',@Yinc_Callback);
        H.firstthick = uicontrol('Parent',uipanel_model,'style','edit','string','100','backgroundcolor','white','Position',[127 159 30 22],'Callback',@firstthick_Callback);
        H.nl = uicontrol('Parent',uipanel_model,'style','edit','string','33','backgroundcolor','white','Position',[127 135 30 22],'Callback',@nl_Callback);
        H.Zinc = uicontrol('Parent',uipanel_model,'style','edit','string','1.3','backgroundcolor','white','Position',[127 111 30 22],'Callback',@Zinc_Callback);
        H.airspacing = uicontrol('Parent',uipanel_model,'style','edit','string','20','backgroundcolor','white','Position',[127 87 30 22],'Callback',@airspacing_Callback);
        
        uicontrol('Parent',uipanel_model,'style','pushbutton','string','Generate Model','Position',[5 46 88 22],'Callback',@M3_generate_model);
        uicontrol('Parent',uipanel_model,'style','pushbutton','string','View Model','Position',[98 46 70 22],'Callback',@M3_view_model);
        uicontrol('Parent',uipanel_model,'style','pushbutton','string','Save as ModEM','Position',[5 24 88 22],'Callback',@M3_save_modem_mod);
        uicontrol('Parent',uipanel_model,'style','pushbutton','string','Save as WS','Position',[98 24 70 22],'Callback',@M3_save_ws_mod);
        uicontrol('Parent',uipanel_model,'style','pushbutton','string','Add Topo','Position',[98 2 70 22],'Callback',@M3_add_topo);
        
    uipanel_edit_model = uipanel('Title','Edit Model','FontSize',8,'units','pixels','Position',[928 205 178 144]);
    % This panel involves editing the model.  
    
        uicontrol('Parent',uipanel_edit_model,'style','pushbutton','string','Add Row','Position',[6 95 78 22],'Callback',@M3_add_row)
        uicontrol('Parent',uipanel_edit_model,'style','pushbutton','string','Add Column','Position',[90 95 78 22],'Callback',@M3_add_col)
        uicontrol('Parent',uipanel_edit_model,'style','pushbutton','string','Delete Row','Position',[6 71 78 22],'Callback',@M3_del_row)
        uicontrol('Parent',uipanel_edit_model,'style','pushbutton','string','Delete Column','Position',[90 71 78 22],'Callback',@M3_del_col)
%         uicontrol('Parent',uipanel_edit_model,'style','pushbutton','string','Add Layer','Position',[6 47 78 22],'Callback',@add_lay_Callback)
%         uicontrol('Parent',uipanel_edit_model,'style','pushbutton','string','Delete Layer','Position',[90 47 78 22],'Callback',@del_lay_Callback)
        uicontrol('Parent',uipanel_edit_model,'style','pushbutton','string','Prev. Layer','Position',[6 23 78 22],'Callback',@M3_prev_layer)
        uicontrol('Parent',uipanel_edit_model,'style','pushbutton','string','Next Layer','Position',[90 23 78 22],'Callback',@M3_next_layer)
        uicontrol('Parent',uipanel_edit_model,'style','pushbutton','string','Jump to Layer','Position',[6 -1 78 22],'Callback',@M3_jump_to_layer)
        H.lay_num = uicontrol('Parent',uipanel_edit_model,'style','edit','string','1','backgroundcolor','white','Position',[90 -1 30 22],'Callback',@lay_num_Callback);
        
        
    uipanel_outputdata = uipanel('Title','Output Data','FontSize',8,'units','pixels','Position',[734 240 190 194]);
    % This panel involves selecting data periods, error floors, and saving data files.
    
        H.min_per = uicontrol('Parent',uipanel_outputdata,'style','edit','string','0.001','backgroundcolor','white','tooltip',sprintf('Enter minimum period to be exported.\nNo periods below this value will be exported.'),'Position',[6 137 36 22],'Callback',@min_per_Callback);
        H.max_per = uicontrol('Parent',uipanel_outputdata,'style','edit','string','1000','backgroundcolor','white','tooltip',sprintf('Enter maximum period to be exported.\nNo periods above this value will be exported.'),'Position',[66 137 36 22],'Callback',@max_per_Callback);
        H.per_skip = uicontrol('Parent',uipanel_outputdata,'style','edit','string','1','backgroundcolor','white','tooltip',sprintf('Export every nth period.\nExample: "1" takes all periods, "2" takes every other...'),'Position',[46 137 18 22],'Callback',@per_skip_Callback);
        H.errflr_Zdiag = uicontrol('Parent',uipanel_outputdata,'style','edit','string','10','backgroundcolor','white','Position',[112 137 18 22],'Callback',@errflr_Zdiag_Callback);
        H.errflr_Zodiag = uicontrol('Parent',uipanel_outputdata,'style','edit','string','10','backgroundcolor','white','Position',[112 113 18 22],'Callback',@errflr_Zodiag_Callback);
        H.errflr_type = uicontrol('Parent',uipanel_outputdata,'style','edit','string','1','backgroundcolor','white','tooltip',sprintf('Type 1: All Z components error = %% of sqrt(|Zxy*Zyx|)\nType 2: Zxx error = %% of |Zxy|, and so on'),'Position',[112 89 18 22],'Callback',@errflr_type_Callback);
        H.errflr_tip = uicontrol('Parent',uipanel_outputdata,'style','edit','string','0.02','backgroundcolor','white','Position',[110 48 28 22],'Callback',@errflr_tip_Callback);
        
        uicontrol('Parent',uipanel_outputdata,'style','text','string','Periods','Position',[30 160 50 15]);
        uicontrol('Parent',uipanel_outputdata,'style','text','string','% Error Floor','Position',[110 160 70 15]);
        uicontrol('Parent',uipanel_outputdata,'style','text','string','Zxx Zyy','Position',[132 132 43 23]);
        uicontrol('Parent',uipanel_outputdata,'style','text','string','Zxy Zyx','Position',[132 108 43 23]);
        uicontrol('Parent',uipanel_outputdata,'style','text','string','Err flr type','tooltip',sprintf('Type 1: All Z components error = %% of sqrt(|Zxy*Zyx|)\nType 2: Zxx error = %% of |Zxy|, and so on'),'Position',[130 84 60 23]);
        uicontrol('Parent',uipanel_outputdata,'style','text','string','Tip err','Position',[138 44 42 23]);
        
        H.resp = uicontrol('Parent',uipanel_outputdata,'style','popupmenu','string',{'Full Tensor (8)';'FT + Tip (12)';'Off-Diag Z Only (4)';'Tipper Only (4)'},'Position',[4 0 95 30],'Callback',@num_resp_Callback); 
        H.conj_Z = uicontrol('Parent',uipanel_outputdata,'style','checkbox','string','Conj. Z','tooltip',sprintf('WARNING: matfile time dep. is e^+iwt,\n wsinv3dmt needs e^-iwt, and ModEM can use either one.'),'Position',[106 22 60 23],'Callback',@conj_Z_Callback);
        H.conj_tip = uicontrol('Parent',uipanel_outputdata,'style','checkbox','string','Conj. Tip','Position',[106 2 60 23],'Callback',@conj_tip_Callback);
        
        uicontrol('Parent',uipanel_outputdata,'style','pushbutton','string','Save WS3D Data','Position',[6 88 92 22],'Callback',@M3_save_ws_data);
        uicontrol('Parent',uipanel_outputdata,'style','pushbutton','string','Save ModEM Data','Position',[6 65 92 22],'Callback',@M3_save_modem_data);
        uicontrol('Parent',uipanel_outputdata,'style','pushbutton','string','View Data','Position',[6 42 66 22],'Callback',@M3_view_data);
        
    uipanel_ocean_model = uipanel('Title','Ocean Model','FontSize',8,'units','pixels','Position',[596 667 128 89]);
        H.ocean_res = uicontrol('Parent',uipanel_ocean_model,'style','edit','string','0.3','backgroundcolor','white','Position',[10 26 22 20],'Callback',@ocean_res_Callback);
        H.fix_ocean = uicontrol('Parent',uipanel_ocean_model,'style','checkbox','string','Fix Ocean','Position',[10 6 68 17],'Callback',@fix_ocean_Callback);
        
        uicontrol('Parent',uipanel_ocean_model,'style','text','string','Ocean Resistivity','Position',[35 28 85 16]);
        
    H.model_text = uicontrol('style','text','string','Mesh info','FontSize',10,'Position',[120 685 280 80]);
    H.data_text = uicontrol('style','text','string','Number of periods = ','FontSize',10,'Position',[405 729 180 20]);
    
    uicontrol('style','text','string','Resistivity','Position',[863 220 53 15])
    uicontrol('style','text','string','Min','Position',[860 207 28 15])
    uicontrol('style','text','string','Max','Position',[892 207 28 15])
    
    H.min_res = uicontrol('style','edit','string','1','backgroundcolor','white','Position',[860 186 28 21],'Callback',@M3_min_res);
    H.max_res = uicontrol('style','edit','string','1000','backgroundcolor','white','Position',[890 186 31 21],'Callback',@M3_max_res);        
    
    
    H.gridlines = uicontrol('style','checkbox','string','No Grid','Position',[834 160 60 17],'Callback',@gridlines_Callback);
    H.station_names = uicontrol('style','checkbox','string','Station Names','Position',[834 140 90 17],'Callback',@station_names_Callback);
       
    uicontrol('style','text','string','Plot Buffer (km)','Position',[910 115 80 21])
    stat_tol = uicontrol('style','edit','string','4','Position',[994 115 17 21],'Callback',@stat_tol_Callback);
    uicontrol('style','text','string','Default Resistivity','Position',[890 90 100 21])
    H.res = uicontrol('style','edit','string','100','Position',[994 90 31 21],'Callback',@res_Callback);          
    uicontrol('style','text','string','Rotation Angle','tooltip','Enter strike angle, positive = clockwise from north','Position',[890 60 100 21])
    H.rotang = uicontrol('style','edit','string','0','tooltip','Enter strike angle, positive = clockwise from north','Position',[994 60 31 21],'Callback',@rotang_Callback);
    uicontrol('style','pushbutton','string','ROTATE STATIONS AND DATA','tooltip',sprintf('Enter angle, positive = clockwise from north\nData are rotated to this angle\nMesh x-direction (North) will be aligned to this angle'),'Position',[900 36 164 22],'Callback',@M3_set_rot);
    uicontrol('style','pushbutton','string','Smooth Last Air Layer','tooltip',(''),'Position',[940 158 164 22],'Callback',@M3_smooth_air);
    uicontrol('style','pushbutton','string','Delete Layers','tooltip',(''),'Position',[940 180 164 22],'Callback',@M3_delete_layers);
    
    % Update H structure, provide default values
    H.axes1 = axes1;
    H.stat_tol = get(stat_tol,'string');
    
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
    if isfield(H,'XX')
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
    if isfield(H,'x')
        if get(hObject,'Value')
            H.tx = text(H.x,H.y,ones(size(H.x)),H.d.site,'fontsize',11,'FontWeight','bold');
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
    clear all; close all; M3D_v1_old_panels
    delete undo*
end

function exit_Callback(~,~)
    close all; clear all
    delete undo*
end

end % end M3