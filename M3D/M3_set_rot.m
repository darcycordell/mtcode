function M3_set_rot(hObject, ~, ~)
    H=guidata(hObject);
    
    if ~isfield(H,'d')
        warndlg('Load data with stations first')
        return
    end
    
    rot = str2double(get(H.rotang,'string')); % angle (in degrees) positive to rotate CW (negative is CCW)
    
    rotchoice = menu(['Station Data will be rotated to ',num2str(rot),' degrees and mesh x-direction will align to ',num2str(rot),' degrees'],'Continue','Cancel');
    
    if rotchoice ==1
        midx = min(H.x_orig) + abs(max(H.x_orig)-min(H.x_orig))/2; % in m
        midy = min(H.y_orig) + abs(max(H.y_orig)-min(H.y_orig))/2;
        
        c = cosd(-rot);    s = sind(-rot); % stations get rotated opposite direction of the data
        rotmat=[c s;-s c];
        rotated=rotmat*[ H.y_orig'-midy; H.x_orig'-midx];
        H.d.y = (rotated(1,:)+midy)'; H.d.x = (rotated(2,:)+midx)';
        
        ry = min(H.d.y) + abs(max(H.d.y)-min(H.d.y))/2; % find amount shifted from origin
        rx = min(H.d.x) + abs(max(H.d.x)-min(H.d.x))/2; 
        
        H.d.y = H.d.y - ry; % get the stations back to centered at x,y = 0
        H.d.x = H.d.x - rx;

        [H.d]=rotate_d(H.d,rot); % need to rotate data the opposite direction of the stations so that x-direction lines up with mesh
                                 % this function rotates data only
                                 
         H.d.Zerr = H.d.Zerr.*0+(1e-9 + 1i*1e-9); % make errors very small to force error floor to be applied
         H.d.tiperr = H.d.tiperr.*0+(1e-9 + 1i*1e-9);
         H.d.rhoerr = H.d.rhoerr.*0 + 1e-9;
         H.d.phaerr = H.d.phaerr.*0 + 1e-9;                                                                 
                                 
        if isfield(H,'m')
            M3_update_model_plot(H)
            warndlg('The current mesh may no longer be compatible with the rotated stations!','Warning')
        else        
            sttol=str2double(H.stat_tol);
            plot(H.d.y./1000,H.d.x./1000,'vk','markerfacecolor','b'); axis equal; xlabel('Y (km)'); ylabel('X (km)')
            axis([min(H.d.y./1000)-sttol max(H.d.y./1000)+sttol min(H.d.x./1000)-sttol max(H.d.x./1000)+sttol]);
            text(H.d.y./1000,H.d.x./1000,H.d.site,'fontsize',11,'FontWeight','bold');
        end

        helpdlg(['Mesh x-direction rotated to ' num2str(rot),' degrees and data rotated to ',num2str(rot),' degrees. Data errors reset! Use error floors.'],'Stations Rotated!')   
    end
    
    H.mesh_rot = rot;
    guidata(hObject, H);
end