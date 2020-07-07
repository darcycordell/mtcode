function plot_rho_z_model(m,d)
% Function which plots the resistivity versus depth curve beneath a chosen
% model cell. User chooses model cell by clicking.
%
% Usage: plot_rho_z_model(m,d) OR plot_rho_z_model(m)
%
% "m" is the model structure
% "d" is OPTIONAL data structure for plotting site locations


u = user_defaults;

%If no data is supplied, then make a NaN data
if ~exist('d','var')
    d = make_nan_data;
end

zax = [u.zmin u.zmax];
col={'b','m','c','g','y','k'}; g = 1; is =1;

close all
set_figure_size(1);
h(1)=subplot(1,2,1);
plot_slice(m,is,d); hold on;
title(['Depth = ',num2str(m.cz(is)/1000),' km b.s.l.']);

while 1
    
    main_menu = menu('','Choose Slice','Choose Points','Clear Figure','Z Limits (km)','Rho Limits (Ohm m)','Add Color Group','Quit');

    
    if main_menu == 1  %Choose depth slice to plot---------------------------
        %Choose depth slice to plot
        prompt={'Depth to plot slice (km)'}; 
        dlg_title='Plot Vertical Slice';
        def={num2str(m.cz(1)/1000)}; %default is the midpoint depth
        num_lines=1;
        dinp = inputdlg(prompt,dlg_title,num_lines,def);
        dinp=str2double(dinp{1});
        tmp = abs(m.cz-dinp*1000);
        [~, is] = min(tmp); %model cell index of closest depth

        h(1)=subplot(1,2,1); cla(h)
        plot_slice(m,is,d); hold on;
        title(['Depth = ',num2str(m.cz(is)/1000),' km b.s.l.']);

    elseif main_menu == 2 %Select points to plot curve---------------------
        
        click=1;disp('Right click to stop');
        while click ~=3
            [hy,hx,click]=ginput(1); %select polygons on horizontal slice
            if click == 3
                break
            end

            h(1)=subplot(1,2,1);
            plot(hy,hx,[col{g},'o'])

            disp(['Points clicked .... x points: [', num2str(hx),'], y points: [',num2str(hy),']']);

            tmp = abs(m.cx/1000-hx);
            [~, indx] = min(tmp); %model cell index of closest depth
            tmp = abs(m.cy/1000-hy);
            [~, indy] = min(tmp); %model cell index of closest depth

            h(2)=subplot(1,2,2); semilogx(squeeze(m.A(indx,indy,:)), m.cz/1000,[col{g},'-']); hold on; axis ij;
            axis([10^u.colim(1) 10^u.colim(2) zax(1) zax(2)])
            axis ij; xlabel('Resistivity (Ohm m)');ylabel('Depth (km)');
            title(['Resistivity Sounding Beneath Point [',num2str(m.cy(indy)),', ',num2str(m.cx(indx)),']']);
        end

    elseif main_menu==3 %Clear figure--------------------------------------
        cla(h(1));cla(h(2));figure(1);h(2)=subplot(1,2,2);title('');
        figure(1);h(1)=subplot(1,2,1);
        plot_slice(m,is,d); hold on;
        title(['Depth = ',num2str(m.cz(is)/1000),' km b.s.l.']);
        g=1;
           
    elseif main_menu==4    %Adjust Z limits -------------------------------
        prompt = {'Minimum depth in kms:','Maximum depth in kms:'};
        dlg_title = 'Change depth axis';
        num_lines = 1;
        def = {num2str(zax(1)),num2str(zax(2))};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        zax=[str2double(answer{1}) str2double(answer{2})]; 
        figure(1);h(2)=subplot(1,2,2);
        axis([10^u.colim(1) 10^u.colim(2) zax(1) zax(2)]);
        
    elseif main_menu==5 %Rho limits----------------------------------------
        prompt = {'Minimum rho in ohm m:','Maximum rho in ohm m:'};
        dlg_title = 'Change resistivity axis';
        num_lines = 1;
        def = {num2str(10^u.colim(1)),num2str(10^u.colim(2))};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        u.colim=[log10(str2double(answer{1})), log10(str2double(answer{2}))];
        figure(1);h(2)=subplot(1,2,2);
        axis([10^u.colim(1) 10^u.colim(2) zax(1) zax(2)]);
        
    elseif main_menu==6 %Add a color group
        
        disp('Add a different color on the line plot or station plot. Maximum 6 groups')
        g=g+1;
        if g>6
            disp('Error: Only; 6 groups possible')
            g=6;
        end
    else
            close all
            break

    end
end

end