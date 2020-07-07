function plot_isosurface(m,d)
% Function which plots the isosurface of a model volume in 3D view
%
% Usage: plot_isosurface(m) OR plot_isosurface(m,d)
%
% "m" is the model structure
% "d" is an OPTIONAL data structure used only to plot site locations and
% geoboundaries
%

u = user_defaults;
[L] = load_geoboundary_file_list;
irun=1;
irep = 0; 
imenu=menu('Isosurface plot options','Plot surface equal to resistivity value','Plot two surfaces','Return');

%Determine x,y and z limits for plotting
if length(u.xylims) ~= 4
    xind = m.npad(1)+1:m.nx-m.npad(1);
    yind = m.npad(2)+1:m.ny-m.npad(2);
else
    xind = nearestpoint(u.xylims(1)*1000,m.cx):nearestpoint(u.xylims(2)*1000,m.cx);
    yind = nearestpoint(u.xylims(3)*1000,m.cy):nearestpoint(u.xylims(4)*1000,m.cy);
end

zind = 1:nearestpoint(u.zmax*1000,m.cz);
zind_orig = zind;

while irun==1
    dlg_title='';
    num_lines=1;
    %Enter isosurface to plot
    if imenu==1 % plotting less than or equal to value
        prompt={'Plot surface less than or equal to (Ohm-m):'}; 
        if irep == 0
            def_res={'100'}; %default is office to Tim's
        end
        dinp = inputdlg(prompt,dlg_title,num_lines,def_res);
    elseif imenu ==2 % plotting two isosurfaces
        prompt = {'Surface 1: less than or equal to:';'Surface 2: greater than or equal to:'};
        if irep == 0
            def_res = {'10';'100'};
        end
        dinp = inputdlg(prompt,dlg_title,num_lines,def_res);        
        if str2double(dinp{2}) < str2double(dinp{1})
            disp('Resistivity of Surface 2 must be greater than resistivity of Surface 1')
            return
        end
    else
        return
    end     
    
    def_res = dinp; % save default entries in case loop is rerun
    B=m.A(xind,yind,zind); 
    %Create a 3D meshgrid over the indices of interest (no padding)
    [X,Y,Z] = meshgrid(m.cy(yind)/1000,m.cx(xind)/1000,m.cz(zind)/1000);
    set_figure_size(1);
    clf
    if exist('d','var') %If data exists, plot it
        plot3(d.y./1000,d.x./1000,d.z./1000,'k.','markersize',12);hold on
        plot_topo(d,2);
        plot_geoboundaries(L,d.origin,d.z);
    end
    
    if imenu == 1 % plot surface less than or equal to value
        res=str2double(dinp{1});
        handle = patch(isosurface(X,Y,Z,B,res)); %patch isosurface
        set(handle, 'FaceColor', [1 0.2 0]', 'EdgeColor', 'none') % make surface red
        title(['Isosurface of ',num2str(res),' ohm-m']); xlabel('Easting (km)'); ylabel('Northing (km)'); zlabel('Elevation b.s.l. (km)')
    else % plot two surfaces
        res1 = str2double(dinp{1}); res2 = str2double(dinp{2}); 
        handle1 = patch(isosurface(X,Y,Z,B,res1)); %patch isosurface
        set(handle1, 'FaceColor', [1 0.2 0]', 'EdgeColor', 'none') % make surface red
        hold on
        B = -B;
        res2 = -res2;
        handle2 = patch(isosurface(X,Y,Z,B,res2)); %patch isosurface
        set(handle2, 'FaceColor', [0.2 0.2 1]', 'EdgeColor', 'none') % make surface blue
        title(['Isosurfaces of ',num2str(-res2),' ohm-m (blue), and ',num2str(res1),' ohm-m (red)']); xlabel('Easting (km)'); ylabel('Northing (km)'); zlabel('Elevation (km a.s.l.)')
    end
    
    view(3); set(gca,'dataaspectratio',[1 1 1/u.ve])
    set(gca,'Zdir','reverse')
    camlight;  camlight(-80,10); lighting gouraud;
    hold on; 

    grid on; box on; %axis equal; %axis tight; 
    axis([m.cy(yind(1)) m.cy(yind(end)) m.cx(xind(1)) m.cx(xind(end)) m.z(zind_orig(1)) m.z(zind_orig(end))]/1000)
    
    %Main menu to loop
    imenu2=menu('','Repeat with different isosurface contour','Change Isosurface Z limits','Return');

    if imenu2==2 %Change z limits
        prompt = {'Minimum depth in kms:','Maximum depth in kms:'};
        dlg_title = 'Change depth axis';
        num_lines = 1;
        def = {num2str(m.cz(zind(1))/1000),num2str(m.cz(zind(end))/1000)};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        minz=str2double(answer{1});
        maxz=str2double(answer{2});
        zind=find(m.cz<=maxz*1000 & m.cz>=minz*1000);
        irep = 1;
    elseif imenu2==3 %Quit
        irun=0;
    else
        irep = 1;
        close all
    end
   
end


end