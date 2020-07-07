function edit_model_2d(m)

u = user_defaults;

close all
set_figure_size(1); hold on
%Set the ranges to be plotted on the cross-section
yind = nearestpoint(u.xylims(3)*1000,m.cy):nearestpoint(u.xylims(4)*1000,m.cy);
zind = 1:nearestpoint(u.zmax*1000,m.cz);

plot_cross_section(m.cy(yind)/1000,m.cz(zind)/1000,m.A(1,yind,zind));   

me = m;

[Z,Y] = meshgrid(m.cz,m.cy);

main_menu = menu('','Add an irregular polygon by clicking','Manually edit point by point','Save Model');

if main_menu == 1 %IRREGULAR POLYGON
    figure(1);
    [hx,hy]=ginputExtra(100); %select polygons on horizontal slice

    disp(['Points clicked .... EW points: [', num2str(hx),'], NS points: [',num2str(hy),']']);

    if isempty(hx) == 1 || length(hx) < 3 %If only 1 or 2 points are clicked, a polygon cannot be made!
        disp('Error: Polygon must be defined by at least 3 points!')

    else

        %Get the row and column indices of the polygon
        [in, on] = inpolygon(Y/1000,Z/1000,hx,hy);
        [col, row]=find(in==1);

        if all(on) == 0 %If on ~=0, then one of the points clicked is EXACTLY on a mesh vertex.
                        %Very unlikely, but it will cause problems.

            %Get the new resistivity value and the depth interval over
            %which to add the polygon
            prompt={'Replace values inside polygon with rho ='};
            dlg_title='Replace rho';
            def={'30'};
            num_lines=1;
            new_rho = inputdlg(prompt,dlg_title,num_lines,def);

            %Replace values inside polygon over the depth range specified
            vol = 0;
            for i = 1:length(row)
                 vol = (me.dy(col(i))/1000)*(me.dz(row(i))/1000)+vol;
                 me.A(1,col(i),row(i)) = str2double(new_rho{1});

                 ind = sub2ind(size(me.A),1,col(i),row(i));
                 R(ind,:) = [1; str2double(new_rho{1})];

            end

            disp(['Volume replaced = ',num2str(vol),' cubic km'])

            set_figure_size(1); hold on
            plot_cross_section(m.cy(yind)/1000,m.cz(zind)/1000,me.A(1,yind,zind));  
        end

    end
elseif main_menu==2 %MANUAL EDIT-------------------------------------------

    %Manually edit model point-by-point by clicking

    i = 1;

    new_rho = 30;

    %Get new resistivity value
    prompt={'Replace values inside polygon with rho ='};
    dlg_title='Replace rho';
    def={num2str(new_rho)};
    num_lines=1;
    new_rho = str2double(inputdlg(prompt,dlg_title,num_lines,def));

    while 1 %Click points on the map to edit cells. Right click to stop
         figure(1);
         [hx,hy,click]=ginput(1); %select polygons on horizontal slice

        if click == 3
            break
        else

            ix = nearestpoint(hx,me.cy/1000,'previous');
            iy = nearestpoint(hy,me.cz/1000,'previous');

            i = i+1;
            me.A(1,ix,iy) = new_rho; hold off

            ind = sub2ind(size(me.A),1,ix,iy);
            R(ind,:) = [1; new_rho];

            set_figure_size(1); hold on
            [~,~,~] = plot_cross_section(m.cy(yind)/1000,m.cz(zind)/1000,me.A(1,yind,zind));  
        end
    end


end