function me = edit_model(m,d)
% Function allows user to edit a model by adding blocks, polygons,
% cell-by-cell clicking, or replace entire layers or columns with a
% specified resistivity value
%
% Usage: edit_model(m,d) OR edit_model(m)
%
% "m" is the model structure
% "d" is an OPTIONAL data structure. It is only used for plotting site
% locations on figures
%
% Note that if no data structure is supplied then it is currently not
% possible to plot a geoboundary on the plots because the geoboundary
% requires a lat-long reference
%
% When saving the model, you can output as a WSINV or ModEM model.
% The code also outputs a mat file which contains the indices and
% resistivity values that were replaced. This can then be loaded back in at
% any time if you want to repeat the same edits on a different model.
%
% Whenever replacing values at indices, the user can enter 'x' when asked
% to input the new resistivity value to return those indices to their
% original model values.
%
% (C) 2020 Unsworth Research Group (University of Alberta, Edmonton, Canada)


%%
%debugging purposes
% clear all
% m = load_model_modem('ModEM_NLCG_200.rho');
% u = user_defaults;
% d = load_data_modem('ModEM_NLCG_200.dat');
%---------------------------------------------

close all

if ~exist('d','var')
    [d] = make_nan_data;
end

save_count = 1; %If multiple models are saved and output in the same session, then file names are given a counter

indz = round(m.nz/2); %midpoint in depth

R = zeros(m.nx*m.ny*m.nz,2); %This variable keeps track of every index that is edited in the model
    %This is done in linear indices. The first column is all the linear
    %indices. A 1 in that column means that the cells has been edited, a 0
    %means the cell has not been edited. The second column gives the new
    %rho value for the edited cells. Note that if you do multiple edits,
    %this variable will keep track of the sum of all edits.

set_figure_size(1);
plot_slice(m,indz,d);

me = m; %"me" is the model_to_edit structure. Keep "m" in memory to start over.

%Some necessary defaults
north = 1; south = -1; east = 1; west = -1; top = 0; bottom = 5; % defaults for model coords
north_ll = -20; south_ll = -30; east_ll = -67; west_ll = -68; % defaults for lat/lon
%%
while 1

% Main Menu
main_menu=menu('','View Model','Add a Rectangular Prism with NS, EW, Depth range in km','Add a Rectangular Prism with Lat, Lon, Depth range','Add an irregular polygon by clicking','Manually edit point by point','Replace Sections/Layers','Load Indices to Replace','Start Over','Save Model','Save Edited Indices','Quit');

if main_menu == 1 %VIEW MODEL----------------------------------------------
    
    view_menu = menu('','Plot Slices','Jump to Layer','Cross Section');
    
    if view_menu == 1 %Loops through the model slice by slice
        set_figure_size(1);
        [indz] = plot_slice_menus(me,indz,d);
    elseif view_menu == 2 % Jump to layer
        prompt={'Layer Depth (km):'};
        dlg_title='Jump to Layer';
        def={num2str(m.z(indz)/1000)};
        num_lines=1;
        dinp = inputdlg(prompt,dlg_title,num_lines,def);

        if ~isempty(dinp)
            indz = nearestpoint(str2double(dinp{1})*1000,m.z);
        end

        set_figure_size(1);
        plot_slice(me,indz,d);
    elseif view_menu == 3 %Cross section
        
        plot_diagonal_section(me,d);
               
    end
        
elseif main_menu == 2 %ADD RECTANGULAR PRISM WITH NS EW DEPTH RANGES IN KM-
    
    %Pick slices manually by entering x-y distances. This is useful for repeatability. Once
    %you have clicked around your conductors in the first option, you can
    %enter the slices and recreate the same model
    prompt={'North Edge (km)','South Edge (km)','East Edge (km)','West Edge (km)','Top (km)','Bottom (km)'};
    dlg_title='Edit Block with edges at N, S, E, W, Top, Bottom';
    def={num2str(north),num2str(south),num2str(east),num2str(west),num2str(top),num2str(bottom)};
    num_lines=1;
    dinp = inputdlg(prompt,dlg_title,num_lines,def);
    
    if ~isempty(dinp)

        %Pull out min and max ranges entered
        north=str2double(dinp{1});
        south=str2double(dinp{2});
        east=str2double(dinp{3});
        west=str2double(dinp{4});
        top=str2double(dinp{5});
        bottom=str2double(dinp{6});

        if north<=south || east<=west || bottom<=top %Check that user inputs are ok
            disp('One of your ranges is incorrect')
            error('When entering NS, EW and top/bottom ranges, enter in model coordinates. North > south, east > west, bottom > top.')
        else
            
            %Get max and min indices of the block to add
            iminy = nearestpoint(west*1000,me.y,'next');
            imaxy = nearestpoint(east*1000,me.y,'previous');
            iminx = nearestpoint(south*1000,me.x,'next');
            imaxx = nearestpoint(north*1000,me.x,'previous');
            iminz = nearestpoint(top*1000,me.z,'next');
            imaxz = nearestpoint(bottom*1000,me.z,'previous')-1;

            if imaxx > m.nx
                imaxx = m.nx;
            end

            if imaxy > m.ny
                imaxy = m.ny;
            end

            %Edit either conductors, resistors or replace box directly
            [me.A,R]=edit_block(R,iminx,imaxx,iminy,imaxy,iminz,imaxz,me.A);

            set_figure_size(1);
            plot_slice(me,iminz,d);
            
            indz = iminz;

        end
    
    end
    
elseif main_menu == 3 %ADD RECTANGULAR PRISM WITH LAT LON DEPTH------------
    
    if ~unique(d.origin) % if no d structure loaded, d.origin = [0 0] from make_nan_data
        disp('No data structure loaded, so latitude and longitude of model are unknown! Load a data structure or edit in model coordinates')
        continue     
    end
    %Pick slices manually by entering lat-lon ranges. Lat-lon are converted
    %to x-y using d.origin
    prompt={'Lat. max (decimal degrees)','Lat. min (decimal degrees)','Lon. max (decimal degrees)','Lon. min (decimal degrees)','Top (km)','Bottom (km)'};
    dlg_title='Edit Block with latitude and longitude ranges (decimal degrees), Top and Bottom in km';
    def={num2str(north_ll),num2str(south_ll),num2str(east_ll),num2str(west_ll),num2str(top),num2str(bottom)};
    num_lines=1;
    dinp = inputdlg(prompt,dlg_title,num_lines,def);
    
    if ~isempty(dinp)

        %Pull out min and max ranges entered
        north=str2double(dinp{1});
        south=str2double(dinp{2});
        east=str2double(dinp{3});
        west=str2double(dinp{4});
        top=str2double(dinp{5});
        bottom=str2double(dinp{6});

        if north<=south || east<=west || bottom<=top %Check that user inputs are ok
            disp('One of your ranges is incorrect')
            error('When entering NS, EW and top/bottom ranges, enter in model coordinates. North > south, east > west, bottom > top.')
        else
                    
            lon = [west east east west];
            lat = [north north south south];
            [y,x] = geo2utm(lon,lat,d.origin(2),d.origin(1)); % follow model convention of y and x
            y = y - 500000;
            west = mean([y(1) y(4)]); % keep edges parallel by taking the mean of the 2 vertices on each edge
            east = mean([y(2) y(3)]);
            south = mean([x(3) x(4)]);
            north = mean([x(1) x(2)]); % these are in m
            
            
            %Get max and min indices of the block to add
            iminy = nearestpoint(west,me.y,'next');
            imaxy = nearestpoint(east,me.y,'previous');
            iminx = nearestpoint(south,me.x,'next');
            imaxx = nearestpoint(north,me.x,'previous');
            iminz = nearestpoint(top*1000,me.z,'next'); % this is still km hence the 1000 factor
            imaxz = nearestpoint(bottom*1000,me.z,'previous')-1;

            %Edit either conductors, resistors or replace box directly
            [me.A,R]=edit_block(R,iminx,imaxx,iminy,imaxy,iminz,imaxz,me.A);

            set_figure_size(1);
            plot_slice(me,iminz,d);
            
            indz = iminz;

        end
    
    end
    
    
elseif main_menu == 4 %ADD POLYGON-----------------------------------------
    
    %Add a polygon by clicking on a given slice and specifying the z range that you want to add the polygon   
    
    %Choose depth slice to plot
    prompt={'Depth to plot slice (km)'}; 
    dlg_title='Plot Vertical Slice';
    def={num2str(me.z(indz)/1000)}; %default is the last ploted slice
    num_lines=1;
    dinp = inputdlg(prompt,dlg_title,num_lines,def);
    
    if ~isempty(dinp)
    
        dinp=str2double(dinp{1});

        indz = nearestpoint(dinp*1000,me.z); %Get index of slice to plot

        set_figure_size(1);
        plot_slice(me,indz,d);

        disp('Click outline of area to paint on horizontal slice. Right click to stop selection.')

        figure(1);
        [hx,hy]=ginputExtra(100); %select polygons on horizontal slice

        disp(['Points clicked .... EW points: [', num2str(hx),'], NS points: [',num2str(hy),']']);

        if isempty(hx) == 1 || length(hx) < 3 %If only 1 or 2 points are clicked, a polygon cannot be made!
            disp('Error: Polygon must be defined by at least 3 points!')

        else
            
            %Get the row and column indices of the polygon
            [in, on] = inpolygon(m.Xc/1000,m.Yc/1000,hx,hy);
            [col, row]=find(in==1);

            if ~all(on) %If on ~=0, then one of the points clicked is EXACTLY on a mesh vertex.
                            %Very unlikely, but it will cause problems.

                %Get the new resistivity value and the depth interval over
                %which to add the polygon
                prompt={'New rho =','Replace value directly (1) OR replace condcutors <new_rho (2)','From depth (km):','To depth (km):'};
                dlg_title='Replace rho';
                def={'30','1',num2str(me.z(indz)/1000),num2str(me.z(indz+1)/1000)};
                num_lines=1;
                new_rho = inputdlg(prompt,dlg_title,num_lines,def);
                
                if ~isempty(new_rho)
                    
                    replace_flag = str2double(new_rho{2});
                    zup = str2double(new_rho{3})*1000;
                    zdown = str2double(new_rho{4})*1000;
                    

                    disp(['Replace rho values with ', new_rho{1},' Ohm m from ', new_rho{2},' km to ',new_rho{3},' km']);

                    %Get z indices
                    isup = nearestpoint(zup,me.z);
                    isdown = nearestpoint(zdown,me.z)-1;

                    %Replace values inside polygon over the depth range specified
                    vol = 0;
                    for i = 1:length(row)
                        for j = isup:isdown
                            if replace_flag==2
                                if ~strcmp(new_rho{1},'x')
                                    if me.A(col(i),row(i),j)<str2double(new_rho{1})
                                        me.A(col(i),row(i),j) = str2double(new_rho{1});
                                        vol = (me.dx(col(i))/1000)*(me.dy(row(i))/1000)*(me.dz(j)/1000)+vol;
                                        
                                        ind = sub2ind(size(me.A),col(i),row(i),j);
                                        R(ind,:) = [1; str2double(new_rho{1})];
                                    end
                                end
                            else
                                 vol = (me.dx(col(i))/1000)*(me.dy(row(i))/1000)*(me.dz(j)/1000)+vol;

                                 if strcmp(new_rho{1},'x') %Special option: if user enters 'x', then the cells inside the polygon will go back to the values of the UNEDITED model
                                    me.A(col(i),row(i),j) =  m.A(col(i),row(i),j); 
                                 else
                                    me.A(col(i),row(i),j) = str2double(new_rho{1});
                                 end

                                 ind = sub2ind(size(me.A),col(i),row(i),j);
                                 R(ind,:) = [1; str2double(new_rho{1})];
                            end

                        end
                    end

                    disp(['Volume replaced = ',num2str(vol),' cubic km'])

                    set_figure_size(1);
                    plot_slice(me,indz,d);
                
                end
            end

        end
    
    end

elseif main_menu==5 %MANUAL EDIT-------------------------------------------
    
    %Manually edit model point-by-point by clicking
    
    i = 1;

    %Plot a slice at a given depth to begin clicking
    prompt={'Depth to plot slice (km)'}; 
    dlg_title='Plot Vertical Slice';
    def={num2str(me.z(indz)/1000)}; %default is the current slice
    num_lines=1;
    dinp = inputdlg(prompt,dlg_title,num_lines,def);
    
    if ~isempty(dinp)
        dinp=str2double(dinp{1});

        %Get z index to plot
        indz = nearestpoint(dinp*1000,me.z);

        set_figure_size(1); hold on
        plot_slice(me,indz,d);

        new_rho = 30;
        while 1 %While loop continues until editing is done. User can select multiple slices

            %Get new resistivity value
            prompt={'Replace values inside polygon with rho ='};
            dlg_title='Replace rho';
            def={num2str(new_rho)};
            num_lines=1;
            inp = inputdlg(prompt,dlg_title,num_lines,def);
            
            if ~isempty(inp)
                if ~strcmp(inp{1},'x')
                    new_rho = str2double(inp{1});
                end

                while 1 %Click points on the map to edit cells. Right click to stop
                     figure(1);
                     [hy,hx,click]=ginput(1); %select polygons on horizontal slice

                    if click == 3
                        break
                    else

                        ix = nearestpoint(hx,me.cx/1000);
                        iy = nearestpoint(hy,me.cy/1000);

                        i = i+1;
                        if strcmp(inp{1},'x') %Special option: if user enters 'x', then the cell clicked will go back to the value of the UNEDITED model
                            new_rho = m.A(ix,iy,indz); hold off
                        end

                        me.A(ix,iy,indz) = new_rho; hold off

                        ind = sub2ind(size(me.A),ix,iy,indz);
                        R(ind,:) = [1; new_rho];

                        clf
                        plot_slice(me,indz,d);
                    end
                end


                %Continue editing on the next or previous slice.
                pmenu = menu('','Next Slice','Previous Slice','Quit');

                if pmenu == 1
                    indz = indz+1;
                    clf
                    plot_slice(me,indz,d);
                elseif pmenu == 2
                    indz = indz-1;
                    clf
                    plot_slice(me,indz,d);
                else
                    break
                end
            
            end
        end
    
    end
    
elseif main_menu == 6 %REPLACE SECTIONS OR LAYERS------------------------------
    
    replace_menu = menu('Replace:','NS Section','EW Section','Layer','Halfspace below','Back');
    
    if replace_menu == 1 % Replace a NS section (i.e. a y-column in the A matrix (or a vertical line in map view))

        %Input the distance in the EW direction to replace and the new rho
        %value
        prompt={'Replace NS Section (km)','New Rho Value = '};
        dlg_title='Replace rho';
        def={'0','30'};
        num_lines=1;
        dinp = inputdlg(prompt,dlg_title,num_lines,def);

        %Find nearest y index
        iy = nearestpoint(str2double(dinp{1})*1000,me.cy);
        new_rho = str2double(dinp{2});
        
        %You must find the topography surface so that the air cells are not
        %replaced
        t = me.Z(:,iy); %Get the topography surface along that section
        for i = 1:length(t)
            %Find the nearest z index to the topography surface
            iz = nearestpoint(t(i),me.cz);
            
            %Replace all z points from the topography surface to the bottom
            %of the model for that row.
            me.A(i,iy,iz:end) = new_rho;
            
            ind = sub2ind(size(me.A),i*ones(length(iz:me.nz),1)',iy*ones(length(iz:me.nz),1)',iz:me.nz);
            R(ind,1) = 1;
            R(ind,2) = new_rho;
            
        end
        
        figure(1)
        plot_slice(me,indz,d)
            
        
    elseif replace_menu == 2 %Replace an EW section (i.e. an x-row in the A matrix (or a horizontal line in map view))
        
        %Input the distance in the NS direction to replace and the new rho
        %value
        prompt={'Replace EW Section (km)','New Rho Value = '};
        dlg_title='Replace rho';
        def={'0','30'};
        num_lines=1;
        dinp = inputdlg(prompt,dlg_title,num_lines,def);

        %Find nearest x index
        ix = nearestpoint(str2double(dinp{1})*1000,me.cx);
        new_rho = str2double(dinp{2});
        
        %You must find the topography surface so that the air cells are not
        %replaced
        t = me.Z(ix,:); %Get the topography surface along that section    
        for i = 1:length(t)
            %Find the nearest z index to the topography surface
            iz = nearestpoint(t(i),me.cz);
            
            %Replace all z points from the topography surface to the bottom
            %of the model for that column.
            me.A(ix,i,iz:end) = new_rho;
            
            ind = sub2ind(size(me.A),ix*ones(length(iz:me.nz),1)',i*ones(length(iz:me.nz),1)',iz:me.nz);
            R(ind,1) = 1;
            R(ind,2) = new_rho;
            
        end
        
        figure(1)
        plot_slice(me,indz,d)
        
    elseif replace_menu == 3 %Replace Layer
        
        %Input the layer to replace and the new rho
        %value
        prompt={'Replace Layer at Depth (km b.s.l.)','New Rho Value = '};
        dlg_title='Replace rho';
        def={'5','30'};
        num_lines=1;
        dinp = inputdlg(prompt,dlg_title,num_lines,def);

        %Get z index of layer
        indz = nearestpoint(str2double(dinp{1})*1000,me.cz);
        new_rho = str2double(dinp{2});
        
        %Find the indices corresponding to the topography surface for that
        %layer
        [tx,ty] = find(me.Z<=me.cz(indz));
        
        if length(tx) ~= me.nx*me.ny %If the topography surface intersects the 
                                     %layer then it is necessary to replace one by one
            for i = 1:length(tx) %Loop over all topography indices
                me.A(tx(i),ty(i),indz) = new_rho; 
                
            end
            
            ind = sub2ind(size(me.A),tx,ty,indz*ones(length(tx),1));
            R(ind,1) = 1;
            R(ind,2) = new_rho;
            
        else %If the topography surface does not intersect the layer, then the replacement
             %can be done quickly with one line
            me.A(:,:,indz) = new_rho; 
            
            r = reshape(R(:,1),size(me.A));
            r(:,:,indz) = 1;
            R(:,1) = r(:);
            
            r = reshape(R(:,2),size(me.A));
            r(:,:,indz) = new_rho;
            R(:,2) = r(:);
            
        end
        
        figure(1);clf
        plot_slice(me,indz,d)
        
        
    elseif replace_menu == 4 %Replace halfspace below

        %Get depth to replace everything below and new rho value
        prompt={'Replace Halfspace Below Depth (km b.s.l.)','New Rho Value = '};
        dlg_title='Replace rho';
        def={'5','30'};
        num_lines=1;
        dinp = inputdlg(prompt,dlg_title,num_lines,def);

        %Get z index of depth to replace
        indz = nearestpoint(str2double(dinp{1})*1000,me.cz);
        new_rho = str2double(dinp{2});
        
        %Find the indices corresponding to the topography surface for that
        %layer
        [t] = find(me.Z<=me.cz(indz));
        
        if length(t) ~= m.nx*m.ny %If the topography surface intersects the 
                                  %layer then it is necessary to replace one by one
            disp('Note: The depth you chose is partially (or wholly) above the topography surface.')
            
            %Find z index of lowest topography point
            indt_max = nearestpoint(max(me.Z(:)),me.cz);
            
            %Replace all points from that z index to the bottom of the
            %model
            me.A(:,:,indt_max:me.nz) = new_rho;

            for i = indz:indt_max-1 %Loop over all the z indices that are above or partially above the topography surface
                [tx,ty] = find(me.Z<=me.cz(i));
                for j = 1:length(tx) %Loop over all topography indices for that layer
                    me.A(tx(j),ty(j),i) = new_rho; 
                end
                
                ind = sub2ind(size(me.A),tx,ty,i*ones(length(tx),1));
                R(ind,1) = 1;
                R(ind,2) = new_rho;
                
            end
            
            
            
        else %If the topography surface does not intersect the layer, then the replacement
             %can be done quickly with one line
            me.A(:,:,indz:me.nz) = new_rho;
            
            r = reshape(R(:,1),size(me.A));
            r(:,:,indz:me.nz) = 1;
            R(:,1) = r(:);
            
            r = reshape(R(:,2),size(me.A));
            r(:,:,indz:me.nz) = new_rho;
            R(:,2) = r(:);
        
        end
        
        figure(1); clf
        plot_slice(me,indz,d)
        
        
    end
    
elseif main_menu == 7 %LOAD INDICES AND REPLACE----------------------------
    %Anytime you save a model, it also outputs the x,y,z,rho indices that
    %were replaced. For repeatability, you can load those indices into the
    %a different model and do the same replacements.
    %
    %This is also useful if you have a different geophysical model (e.g.
    %gravity, seismic) which you have interpolated onto your MT model grid
    %and then grabbed those indices

    [index_filename, filepath]=uigetfile({'*.mat'},'Choose Initial ModEM model file'); 
    curdir = pwd;
    
    if index_filename ~= 0
    %Option here to load those indices and edit
        cd(filepath)
        load(index_filename); cd(curdir);
        
                    
        if exist('ix','var')
            indices(:,1) = ix;
            indices(:,2) = iy;
            indices(:,3) = iz;
            indices(:,4) = 1;
            clear ix iy iz
        end
        
        if length(indices(1,:))~=4
            %Get depth to replace everything below and new rho value
            prompt={'New Rho Value for Indices = '};
            dlg_title='Replace rho';
            def={'30'};
            num_lines=1;
            dinp = inputdlg(prompt,dlg_title,num_lines,def);

            new_rho = str2double(dinp{1});
            
            indices(:,4) = new_rho;
        end
            
        vol = 0;
        for i = 1:length(indices(:,1))
           me.A(indices(i,1),indices(i,2),indices(i,3)) = indices(i,4);
           vol = (me.dx(indices(i,1))/1000)*(me.dy(indices(i,2))/1000)*(me.dz(indices(i,3))/1000)+vol;
        end
        
        figure(1); clf
        plot_slice(me,min(indices(:,3)),d)
        
        ind = sub2ind(size(m.A),indices(:,1),indices(:,2),indices(:,3));
        R(ind,1) = 1;
        R(ind,2) = indices(:,4);
    
    end
    
elseif main_menu == 8 %START OVER------------------------------------------
    %%
    %Resets the model back to the original model
    
    me = m;
    
    R = zeros(m.nx*m.ny*m.nz,2); 
    
    set_figure_size(1);
    plot_slice(m,indz,d);
    
elseif main_menu == 9 %SAVE MODEL------------------------------------------
    %Save model as WSINV or ModEM format
    curdir = pwd;
    outputfile = ['edit_',num2str(save_count),'_',me.name];
    
    [outputfile, filepath]=uiputfile({'*.rho'},'Save file as',outputfile); 
    
    if outputfile ~= 0

        save_menu=menu('','Save as ModEM','Save as WSINV');
        
        [ixr,iyr,izr] = ind2sub(size(me.A),find(R(:,1)==1));
        
        indices = [ixr iyr izr R(R(:,2)~=0,2)];
        
        save([outputfile,'_indices.mat'],'indices');
        
        cd(filepath)
        if save_menu==1
            write_model_modem(outputfile,me.dx,me.dy,me.dz,me.A,me.origin); 
            save_count = save_count+1;
        elseif save_menu==2 %Rewrote save_model function and broke into two separate functions
            write_model_wsinv(['wsinv_',outputfile],me.nx,me.ny,me.nz,me.dx,me.dy,me.dz,me.A,d.niter,0);
            save_count = save_count+1;
        end
        
        cd(curdir)
    
    end
elseif main_menu == 10 %SAVE EDITED INDICES
    
    [ixr,iyr,izr] = ind2sub(size(me.A),find(R(:,1)==1));

    indices = [ixr iyr izr R(R(:,2)~=0,2)];

    save(['indices.mat'],'indices');

else
    close all
    break
end
end

end %END MAIN--------------------------------------------------------------





        
      


    
            
      
