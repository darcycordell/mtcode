function plot_rho_z_station(m,d)
% Function which plots the resistivity versus depth curve beneath a chosen
% MT site. User chooses the sites to plot.
%
% Usage: plot_rho_z_model(m,d)
%
% "m" is the model structure
% "d" is data structure

%------------------Plot a Slice to Choose Diagonal Slice------------------
%%
u = user_defaults;

close all
set_figure_size(1);
h(1) = subplot(1,2,1);
plot_slice(m,1,d); hold on;
title(['Depth = ',num2str(m.cz(1)/1000),' km b.s.l.']);

res = zeros(m.nz,d.ns);
for i = 1:d.ns
    res(:,i)=squeeze(m.A(d.idx(i),d.idy(i),:));
end

zax = [u.zmin u.zmax];
col={'b','m','c','g','y','k'};
is = 1; next = 1;
next_menu=0; g=1;

while next==1
    
    dummy_nxt=next_menu;
    next_menu=menu('','Plot Next Station','Plot All Stations','Select Stations from Map','Select Stations from Curve','Clear Figure','Z Limits (km)','Rho limits (Ohm m)','Add Color Group','Save All Soundings','Quit');              

    if next_menu==1 %Plot next station-------------------------------------
        figure(1);subplot(1,2,1)
        plot(d.y(is)/1000,d.x(is)/1000,'rs','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6);box on;
        %to plot sounding curve
        figure(1);h(2)=subplot(1,2,2);
        semilogx(res(:,is), m.cz/1000,'-r')
        axis([10^u.colim(1) 10^u.colim(2) zax(1) zax(2)])
        axis ij; xlabel('Resistivity (Ohm m)');ylabel('Depth (km)');
        title(['Resistivity Sounding for Site ',d.site{is}]);
        hold on;
        is=is+1;
        if is>d.ns
            is=d.ns;
        elseif is<1
            is=1;
        end
        
    elseif next_menu==2 %Plot all stations---------------------------------
        figure(1);subplot(1,2,1);
        plot(d.y./1000,d.x./1000,'ks','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6); hold on;
        for is2=1:d.ns
            %to plot sounding curve
            figure(1);h(2)=subplot(1,2,2);
            semilogx(res(:,is2), m.cz/1000,'-r')
            axis([10^u.colim(1) 10^u.colim(2) zax(1) zax(2)])
            axis ij; xlabel('Resistivity (Ohm m)');ylabel('Depth (km)');
            title('Resistivity Sounding for Site All Sites');
            hold on
        end
        
    elseif next_menu==3 %Select Station from Map---------------------------
        click=1;disp('Right click to stop');
        while click~=3

            figure(1);subplot(1,2,1);
             %Clickable section:
            [y,x, click] = ginput(1);%get x,y points on map in latitude and longitude
            if click==3
                break
            end
            difx = abs(d.x-x.*1000);dify=abs(d.y-y.*1000);dif_av=(difx+dify)./2;
            [~, id]=min(dif_av);%index of closest point to mouse click
            is = id;
            plot(d.y(is)/1000,d.x(is)/1000,[col{g},'s'],'MarkerEdgeColor',col{g},'MarkerFaceColor',col{g},'MarkerSize',6);box on;
            title(strrep(d.site(is),'_','\_'));
            %to plot sounding curve
            figure(1);h(2)=subplot(1,2,2);
            semilogx(res(:,is), m.cz/1000,['-',col{g}],'LineWidth',1.5)
            axis([10^u.colim(1) 10^u.colim(2) zax(1) zax(2)])
            axis ij; xlabel('Resistivity (Ohm m)');ylabel('Depth (km)');
            title(['Resistivity Sounding for Site ',d.site{is}]);
            hold on;
        end
        
    elseif next_menu==4 %Select station from curve-------------------------
        click=1;disp('Right click to stop');
        while click~=3

            figure(1);subplot(1,2,2);
             %Clickable section:
            [x, y, click] = ginput(1);%get x,y points on map in latitude and longitude
            if click==3
                break
            end
            difx = abs(res-x);dify=abs(m.cz-y*1000);
            dif_dist = zeros(m.nz,d.ns);
            for k=1:d.ns
                dif_dist(:,k)=sqrt(difx(:,k).^2+dify.^2);
            end
            [dif_min]=min(dif_dist);%index of closest point to mouse click
            [~, id]=min(dif_min);
            is = id;
            figure(1);subplot(1,2,1);
            plot(d.y(is)/1000,d.x(is)/1000,[col{g},'s'],'MarkerEdgeColor',col{g},'MarkerFaceColor',col{g},'MarkerSize',6);box on;
            title(strrep(d.site(is),'_','\_'));
            %to plot sounding curve
            figure(1);h(2)=subplot(1,2,2);
            semilogx(res(:,is), m.cz/1000,['-',col{g}],'LineWidth',1.5)
            axis([10^u.colim(1) 10^u.colim(2) zax(1) zax(2)])
            axis ij; xlabel('Resistivity (Ohm m)');ylabel('Depth (km)');
            title(['Resistivity Sounding for Site ',d.site{is}]);
            hold on;
        end


    elseif next_menu==5 %Clear figure--------------------------------------
        if dummy_nxt~=0 && dummy_nxt~=5
            cla(h(2));figure(1);subplot(1,2,2);title('');
            cla(h(1));figure(1);subplot(1,2,1);
            %to plot sites location
            plot_slice(m,1,d); hold on;
            title(['Depth = ',num2str(m.cz(1)/1000),' km b.s.l.']);
            is=1; g=1;
        end    
        
    elseif next_menu==6    %Adjust Z limits -------------------------------
        prompt = {'Minimum depth in kms:','Maximum depth in kms:'};
        dlg_title = 'Change depth axis';
        num_lines = 1;
        def = {num2str(zax(1)),num2str(zax(2))};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        zax=[str2double(answer{1}) str2double(answer{2})]; 
        figure(1);subplot(1,2,2);
        axis([10^u.colim(1) 10^u.colim(2) zax(1) zax(2)]);
        
    elseif next_menu==7 %Rho limits----------------------------------------
        prompt = {'Minimum rho in ohm m:','Maximum rho in ohm m:'};
        dlg_title = 'Change resistivity axis';
        num_lines = 1;
        def = {num2str(10^u.colim(1)),num2str(10^u.colim(2))};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        u.colim=[log10(str2double(answer{1})), log10(str2double(answer{2}))];
        figure(1);subplot(1,2,2);
        axis([10^u.colim(1) 10^u.colim(2) zax(1) zax(2)]);
        
    elseif next_menu==8 %Add Color group-----------------------------------
        disp('Add a different color on the line plot or station plot. Maximum 6 groups')
        g=g+1;
        if g>6
            disp('Error: Only; 6 groups possible')
            g=6;
        end
        
    elseif next_menu == 9 %Save all soundings as a mat file
        
        for i = 1:d.ns
            res(:,i) = m.A(d.idx(i),d.idy(i),:);
        end
        
        depth = m.cz;
        site = d.site;

        save('station_rho_z_curves.mat','res','depth','site')
    else
        next=0;
        close all
    end
end




end