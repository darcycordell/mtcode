function M3_load_modem_mod(hObject, ~, ~)
% Load ModEM format model file

    H=guidata(hObject);
    delete undo*
    
    if ~isfield(H,'dat')
        warndlg('Load data file first')
        return
    end

    set(H.axes1,'HandleVisibility','ON');
    axes(H.axes1);

    [mod_mod, mod_path] = uigetfile('*.*','Select ModEM model file');
    
    [m] = load_model_modem([mod_path mod_mod]);

    H.top = m.origin(3);
    H.nx = m.nx; 
    H.ny = m.ny;
    H.nz = m.nz;
    H.XX = m.dx;
    H.YY = m.dy;
    H.Z = m.dz;
    H.AA = m.A;
    H.lay = 1; % just set view to first layer

    %now set the model variables so they match the generated model variables

    H.AA(H.nx+1,:,:)=H.AA(H.nx,:,:);
    H.AA(:,H.ny+1,:)=H.AA(:,H.ny,:);
    H.AA(:,:,H.nz+1)=H.AA(:,:,H.nz);

    H.XX = (m.x./1000)';
    H.YY = m.y./1000;
    H.Z = (m.z./1000)';
    
    H.nx=H.nx+1;
    H.ny=H.ny+1;
    H.nz=H.nz+1;
    
    H.m = m;

    M3_update_model_plot(H)
        
    guidata(hObject, H);

end % end load_modem_mod_Callback