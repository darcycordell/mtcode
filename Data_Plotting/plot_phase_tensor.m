function plot_phase_tensor(d)
%%
% Function which plots the phase tensors for a specific MT site at 9
% frequencies determined from the data. Also plots the maximum and minimum
% phase tensor axes, alpha, beta skew and geoelectric strike as a function
% of period.
%
% Usage: plot_phase_tensor(d)
%
% "d" is data input in standard format structure

table(d.site);
    
[sel,val] = listdlg('PromptString',['Select stations (',num2str(d.ns),' total) '],'ListString',d.site,'Name','Stations','ListSize',[300 300]);
    
if ~val
    disp('No sites selected, returning to main');
    return
end

is = 1; rot_ang = 0;

while 1
    
    [p] = calc_phase_tensor(d.Z(:,:,sel(is))); %Calculate phase tensor variables and put them in structure "p"
        %This structure contains phimin, phimax, beta, alpha, etc.
    
    set_figure_size(1);
    
    %Determine the 9 periods that will be plotted
    [~, T_first]  = max(~isnan(d.Z(:,2,:)), [], 1); %Find first NaN period index
    T_first = min(T_first);
    
    [~, T_last]  = max(flipud(~isnan(d.Z(:,2,:))), [], 1); %Find last NaN period index
    T_last = d.nf-min(T_last);
    
    index = floor(linspace(T_first,T_last,9)); %Index of periods to plot
    
    n = 1;
    for ifreq = index
        
        ifreq = min(length(d.T),ifreq+2);
        subplot(3,3,n) %A total of 9 subplots will be plotted

        plot(0,0,'.'); hold on;
    
        % Plot whole ellipse   
        plot(p.y(:,ifreq),p.x(:,ifreq),'k-') 
      
        %Normalized radius of ellipse
        r = 1.1*max([squeeze(abs(p.x(:,ifreq))); squeeze(abs(p.y(:,ifreq)))]);

        % Now plot axis showing co-ordinate frame
        xc1 = [0, r*cos(rot_ang*pi/180.)];      yc1 = [0, r*sin(rot_ang*pi/180.)] ;  
        xc2 = [0, r*cos(pi/2+rot_ang*pi/180.)]; yc2 = [0, r*sin(pi/2+rot_ang*pi/180.)];   
        plot(yc1,xc1,'r-'); plot(yc2,xc2,'b-');
        if ~isnan(r)
           axis([-r r -r r]);
           axis equal
           title(['T = ',num2str(d.T(ifreq)),' s'])
        else % if there is no data at this period
           title(['No Data for T = ',num2str(d.T(ifreq)),' s'])
        end

        n=n+1;
    end
    
    annotation('textbox', [0 0.92 1 0.08], ...
    'String', ['Phase Tensors For Site ',d.site{sel(is)}], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center','FontSize',12,'FontWeight','bold')

    print_figure(['phase_tensor_',d.niter],['phase_tensor_site_ellipse',num2str(sel(is),'%03.0f'),'_',d.site{sel(is)}]); %Save figure
    
    %Plot alpha, beta skew, phimin, phimax etc.
    set_figure_size(2);
    subplot(3,1,1)  % phase_min and phase_max
    semilogx(d.T,p.phimin,'r.') ; hold on;
    semilogx(d.T,p.phimax,'b.') ;
    xlabel ('Period (s)'); ylabel ('Phase (Degrees)');
    legend('Phi Min','Phi Max')
    title(['Phase Tensor Parameters for Site  ',d.site{sel(is)}])

    % plot alpha and strike
    subplot(3,1,2)  
    semilogx(d.T,p.alpha,'.r'); hold on
    semilogx(d.T,p.alpha-p.beta,'-k','LineWidth',2);
    xlabel('Period (s)')
    ylabel('Degrees')
    legend('Alpha','Geoelectric Strike (\alpha - \beta )')
    
    % plot beta
    subplot(3,1,3)
    semilogx(d.T,p.beta,'.b'); hold on
    semilogx([d.T(1) d.T(end)],[-3 -3],'--k');
    semilogx([d.T(1) d.T(end)],[3 3],'--k');
    xlabel('Period (s)')
    ylabel('Degrees')
    legend('Beta');
    
    print_figure(['phase_tensor_',d.niter],['phase_tensor_site_',num2str(sel(is),'%03.0f'),'_',d.site{sel(is)}]); %Save figure
    

    next_menu = menu('','Next Station in List','Previous Station in List','Return');

    if next_menu == 1
        is = is+1;
        if is>length(sel)
            is = length(sel);
        end
    elseif next_menu == 2
        is = is-1;
        if is<1
            is = 1;
        end
    else
        close all
        break
    end


end



