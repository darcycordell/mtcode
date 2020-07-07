function plot_conductance(m,d)
%Function calculates the conductance for a given 3D model. It also
%allows the user to specify over which depth range to calculate the
%conductance
%
% Usage: plot_conductance(m,d) OR plot_conductance(m)
%
% "m" is the model structure
% "d" is an OPTIONAL data structure for plotting site locations on the plot
%


u = user_defaults; 
[L] = load_geoboundary_file_list;


close all
crun=1; 
while crun==1 % Exits on quitting
    crun2=1;

    %Find the depth range over which to calculate the conductance
    prompt={'Min depth to calculate conductance from (km):','Max depth to calculate conductance to (km):'}; 
    dlg_title='Conductance Depth';
    def={num2str(min(m.cz)/1000),num2str(max(m.cz)/1000)}; %default is the bottom of the model and top of model
    num_lines=1;
    dinp = inputdlg(prompt,dlg_title,num_lines,def);

    is=[1,1];
    for i = 1:length(dinp)

        [is(i)] = nearestpoint(str2double(dinp{i})*1000,m.cz);
    end

    is = sort(is);
            

    %Calculate the conductance
    A_sig=1./m.A; %conductivity in S/m
    con = zeros(m.nx,m.ny,1);
    for i=1:m.nx
        for j=1:m.ny
            con_cell=zeros(m.nz,1);
            for k=is(1):is(2)-1
                con_cell(k)=A_sig(i,j,k)*(m.cz(k+1)-m.cz(k)); %Conductance of a given cell
            end
            con(i,j,1)=nansum(con_cell); % Conductance beneath a point (i,j)
        end
    end
    
    disp(['Maximum Conductance =',num2str(max(con(:)))])

    mx=median(con(~isnan(con)))+std(con(~isnan(con))); mn=median(con(~isnan(con)))-std(con(~isnan(con))); if mn<0; mn=0; end;
    while crun2==1 %Exits when colorbar is good
        close all
        %Plot Conductance on a map
        m_pcolor(m.LON,m.LAT,con(:,:,1)); shading flat; hold on
        plot_geoboundaries(L);
        
        if exist('d','var')
            m_plot(d.loc(:,2),d.loc(:,1),'k.','markersize',12);
        end
        
        m_grid('box','fancy','xlabeldir','end','tickdir','in');
        colormap(flipud(u.cmap)); caxis([mn mx]); hcb = colorbar;
        hcb.Label.String = 'Conductance (S)';
        title('Conductance Map');
        
        
        %Begin Menu
        cmenu=menu('','Change Colorbar Limits','Change Depth Limits','Quit');
        if cmenu==1 %Change colorbar limits--------------------------------
            prompt={'Min conductance','Max conductance:'}; 
            dlg_title='Conductance Limits';
            def={num2str(mn),num2str(mx)}; %default is the bottom of the model
            num_lines=1;
            dinp = inputdlg(prompt,dlg_title,num_lines,def);
            mx=str2double(dinp{2}); mn=str2double(dinp{1});
            if mn>mx
                disp('Error: Minimum conductance must be less than Maximum conductance')
                crun2=0;
            end
            
        elseif cmenu==2 %Change depth limits over which to compute conductance-------------
            close all
            crun2=0;
        else
            crun2=0; crun=0;
        end
    end
end

   
end