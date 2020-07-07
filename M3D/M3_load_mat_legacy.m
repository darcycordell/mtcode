function M3_load_mat_legacy(hObject,~)
% Load data from .mat file

    H=guidata(gcbo);

    [dat,~] = uigetfile({'*.mat'},'pick mat file');
    if dat == 0
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

    matfile = load(dat);
    H.dat = dat;
    
    d.loc = [matfile.coord(:,2) matfile.coord(:,1) -matfile.coord(:,3)]; % negated z here for depth b.s.l.
    d.Z = permute(matfile.ZZZ,[2 1 3]); 
    d.Zerr = permute(matfile.EEE,[2 1 3]) + (1i .* permute(matfile.EEE,[2 1 3]));
    d.tip = permute(matfile.ZZZ_HZ,[2 1 3]);
    d.tiperr = permute(matfile.EEE_HZ,[2 1 3]) + (1i .* permute(matfile.EEE_HZ,[2 1 3]));
    d.site = matfile.site';
    d.T = matfile.T;
    d.ns = length(d.Z(1,1,:));
    d.nf = length(d.Z(:,1,1));
    d.nr = 4;
    
    H.d = d;
    
    long = H.d.loc(:,2);
    lat  = H.d.loc(:,1); 
    H.cent_lat  = (max(lat)+min(lat))/2;
    H.cent_long = (max(long)+min(long))/2;
    % Convert station coords from degrees to km N-S and E-W
    [H.d.y,H.d.x] = geo2utm(long,lat,H.cent_long,H.cent_lat);
    H.d.y = H.d.y - 500000; % 500 kilometers is subtracted to put the center of the mesh at 0km (UTM uses 500km as the cent_long point)
    
    H.d.z = H.d.loc(:,3); % x = NS direction, y = EW direction, z = depth b.s.l. as in ModEm inversion
    
    H.x_orig=H.d.x; H.y_orig=H.d.y; % not used much, set_rot uses this but maybe not needed
    plot(H.d.y./1000,H.d.x./1000,'vk','markerfacecolor','b'); axis equal; xlabel('Easting (km)'); ylabel('Northing (km)')
    axis([min(H.d.y./1000)-sttol max(H.d.y./1000)+sttol min(H.d.x./1000)-sttol max(H.d.x./1000)+sttol]);
    text(H.d.y./1000,H.d.x./1000,H.d.site,'fontsize',11,'FontWeight','bold');

    guidata(hObject, H);
    
end % end load_mat_Callback