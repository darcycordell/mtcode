function M3_load_modem_data(hObject, ~)
    
    H=guidata(gcbo);
    
    [file,~] = uigetfile({'*.data'},'Pick ModEM data file');
    if file == 0
        return
    end
    
    M3_clear_var_in_H(hObject,H);
    
    H.mesh_rot = 0;

    delete undo*
    H.undo=0;

    H.lay=1;
    sttol=str2double(H.stat_tol);

    set(H.axes1,'HandleVisibility','ON');
    axes(H.axes1)
    
    H.d = load_data_modem(file);
    H.dat = file; % save filename only for having a seed name for the parameter file
    
    % edit 2019 - plot x and y from data file instead. lat and lon in ModEm data file won't be
    % right if file generated from rotated mesh
    long = H.d.loc(:,2);
    lat  = H.d.loc(:,1); 
    H.cent_lat  = (max(lat)+min(lat))/2;
    H.cent_long = (max(long)+min(long))/2;
%     % Convert station coords from degrees to km N-S and E-W
    [H.d.y,H.d.x] = geo2utm(long,lat,H.cent_long,H.cent_lat);
    H.d.y=H.d.y-500000; % 500 kilometers is subtracted to put the center of the mesh at 0km (UTM uses 500km as the cent_long point)

    H.d.z = -H.d.loc(:,3); % note x = NS direction, y = EW direction, z =  depth b.s.l. as in ModEM inversion
%     H.x = H.d.y./1000; H.y = H.d.x./1000; % switch x and y for M3 plotting
    
    H.x_orig=H.d.x; H.y_orig=H.d.y;
    plot(H.d.y./1000,H.d.x./1000,'vk','markerfacecolor','b'); axis equal; xlabel('Easting (km)'); ylabel('Northing (km)')
    axis([min(H.d.y./1000)-sttol max(H.d.y./1000)+sttol min(H.d.x./1000)-sttol max(H.d.x./1000)+sttol]);
    text(H.d.y./1000,H.d.x./1000,H.d.site,'fontsize',11,'FontWeight','bold');
    
    guidata(hObject, H);   
    
end % end load_modem_data_Callback