function M3_update_model_plot(H)
    % Dec 2019 - x and y directions made consistent with ModEM
    min_res = log10(str2double(get(H.min_res,'string'))); % Get the min resistivity
    max_res = log10(str2double(get(H.max_res,'string'))); % Get the max resistivity
    colim=[min_res-0.15 max_res+0.15]; % caxis values for color plots
    %0.15 to make figures same as pmv* because WS did so...
    cmap=jet(24);
    cmap([23 21 19 17 15 13 11 9 7 5 3],:)=[];
    cmap=flipud(cmap);

    set(H.axes1,'HandleVisibility','ON');
    axes(H.axes1);
    if isfield(H,'AAt')
        aa=reshape(H.AAt(:,:,H.lay),H.nx,H.ny);
    else
        aa=reshape(H.AA(:,:,H.lay),H.nx,H.ny);
    end
    aa(aa==1e17)=NaN; % make air cells plot as white with pcolor
    pcolor(H.axes1,H.YY,H.XX,log10(aa)); hold on; 
    plot(H.d.y./1000,H.d.x./1000,'vk','markerfacecolor','b'); hold off
    axis equal; axis tight
    title(['Layer = ',num2str(H.lay),', Depth = ',num2str(H.Z(H.lay)),' km']);
    xlabel('Y (km)'); ylabel('X (km)')
    caxis(colim); colormap(cmap);
    hcb=colorbar('ver','position',[0.7 0.1 0.02309 0.1977]);
    set(hcb,'YTick',[-1,-.5,0 0.5,1,1.5,2,2.5,3,3.5,4],...
        'YTickLabel',{'0.1','0.3','1','3','10','30','100','300','1000','3000','10000'});
    ylabel(hcb,'Resistivity (\Omegam)','fontsize',12);
    set(hcb,'yaxislocation','left');
   
    newstring = sprintf('%s\n%s\n%s\n%s\n%s','Mesh info',['Maximum depth=',num2str(H.Z(end),'%.2f'),' km'],...
        ['E-W Dist from center =',num2str(abs(min(H.YY)),'%.2f'),' km'],...
        ['N-S Dist from center =',num2str(abs(min(H.XX)),'%.2f'),' km'],...
        ['NX=',num2str(H.nx-1),', NY=',num2str(H.ny-1),', NZ=',num2str(H.nz-1)]);
    set(H.model_text,'string',newstring)
    
end