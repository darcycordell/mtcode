function [x,y,z,rho] = plot_cross_section_3d(m,flag,d)
% Function to plot NS or EW cross-sections through a 3D resistivity model
% User chooses the section to plot on the map slice.
%
% Usage: [x,y,z,rho] = plot_cross_section_3d(m,flag,d)
%   OR   [x,y,z,rho] = plot_cross_section_3d(m,flag)
%
% Inputs:
% "m" is the model structure
% "flag" specifies whether the section should be plotted NS or EW
%       flag = 'NS' 
%       flag = 'EW'
%
% "d" is OPTIONAL to data structure to plot site locations on sections
%
% Outputs:
% x,y,z are the NS, EW and depth locations of the section
% rho is the matrix of resistivity values associated with the section

u = user_defaults;
close all
        

%Plot midpoint slice
id = round(m.nz/2);
figure(1);
if exist('d','var')
    plot_slice(m,id,d);
else
    plot_slice(m,id);
end

title(['Depth = ',num2str(m.cz(id)/1000),' km b.s.l.']);

%Get input from user by clicking the slice they would like to plot
figure(1);
[yp,xp]=ginput(1);  %graphically choose the slice to be plotted


%Plot the line of the chosen section on the map view slice
if strcmp(flag,'NS')

    %Plot the chosen NS section on the map
    [idy, ~] = nearestpoint(yp*1000,m.cy);
    plot([m.cy(idy)/1000 m.cy(idy)/1000],[-10^8 10^8],'-m','LineWidth',3)
    
    %Plot the sites that will be projected onto the section with a tolerance
    %distance of u.tol (see user_defaults.m)
    if exist('d','var') %Only if a data structure is supplied
        if isfield(d,'idy')
            site_disty = abs(d.idy-idy);
            plot(d.y(site_disty<=u.tol)/1000,d.x(site_disty<=u.tol)/1000,'m.','markersize',15)
            plot([m.cy(idy)/1000 m.cy(idy)/1000]+u.tol,[-10^8 10^8],'-m')
            plot([m.cy(idy)/1000 m.cy(idy)/1000]-u.tol,[-10^8 10^8],'-m')
        else
           disp('Link geographic and model coordinates to plot station locations on cross-section. See link_model_data.m')
        end
    end


elseif strcmp(flag,'EW')
    
    %Plot the chosen EW section on the map
    [idx, ~] = nearestpoint(xp*1000,m.cx);
    plot([-10^8 10^8],[m.cx(idx)/1000 m.cx(idx)/1000],'-r','LineWidth',3)
    
    %Plot the sites that will be projected onto the section with a tolerance
    %distance of u.tol (see user_defaults.m)
    if exist('d','var') %Only if a data structure is supplied
        if isfield(d,'idx')
            site_distx = abs(d.idx-idx);
            plot(d.y(site_distx<=u.tol)/1000,d.x(site_distx<=u.tol)/1000,'r.','markersize',15)
            plot([-10^8 10^8],[m.cx(idx)/1000 m.cx(idx)/1000]+u.tol,'-r')
            plot([-10^8 10^8],[m.cx(idx)/1000 m.cx(idx)/1000]-u.tol,'-r')
        else
            disp('Link geographic and model coordinates to plot station locations on cross-section. See link_model_data.m')
        end
    end
    
    


end
%
%
%--------------------------------------------------------------------------
%Now Plot the Cross-Section
if strcmp(flag,'NS')
    
    %Set the ranges to be plotted on the cross-section
    
    if length(u.xylims) ~= 4
        xind = m.npad(1)+1:m.nx-m.npad(1);
        axlims = [sort([m.x(min(xind)) m.x(max(xind))])/1000 u.zmin u.zmax];
    else
        xind = nearestpoint(u.xylims(1)*1000,m.cx,'previous'):nearestpoint(u.xylims(2)*1000,m.cx,'next');
        axlims = [u.xylims(1) u.xylims(2) u.zmin u.zmax];
    end
    
    if any(isnan(xind))
        xind = 1:m.nx;
    end
    
    zind = 1:nearestpoint(u.zmax*1000,m.cz,'next');
    
    if any(isnan(zind))
        zind = 1:m.nz;
    end
    
    x = [m.x(xind); m.x(xind(end)+1)]/1000;
    z = [m.z(zind); m.z(zind(end)+1)]/1000;
    C = squeeze(log10(m.A(xind,idy,zind)))';
    C = horzcat(C,C(:,end));
    C = vertcat(C,C(end,:));
    
    figure(2); hold on
    [x,z,rho] = plot_cross_section(x,z,C);
    axis(axlims)
    title(['North-South Section @ ',num2str(m.cy(idy)/1000),' km (EW)']);
    set(gca,'Layer','top')
    
    if strcmp(u.gridlines,'off')
        shading flat
    end
    
    %If a data structure exists, project the nearby sites onto the section
    if exist('d','var') && isfield(d,'idx')
        plot(d.x(site_disty<=u.tol)/1000,m.Z(d.idx(site_disty<=u.tol),idy)/1000 -u.zoff,'vk','MarkerFaceColor','k')
    end
    print_figure(['vertical_profiles_',m.niter],['N-S_at_slice_',num2str(idy,'%03.0f')]);
    
    %Save figure
    figure(1)
    print_figure(['vertical_profiles_',m.niter],['N-S_at_slice_',num2str(idy,'%03.0f'),'_MAP']);
    
    %Output the x,y,z,rho coordinates of the chosen slice
    y = ones(length(x),1).*m.cy(idy)/1000;
    rho = squeeze(rho);

elseif strcmp(flag,'EW')
    
    %Set the ranges to be plotted on the cross-section
    if length(u.xylims) ~= 4
        yind = m.npad(2)+1:m.ny-m.npad(2);
        axlims = [sort([m.y(min(yind)) m.y(max(yind))])/1000 u.zmin u.zmax];
    else
        yind = nearestpoint(u.xylims(3)*1000,m.cy,'previous'):nearestpoint(u.xylims(4)*1000,m.cy,'next');
        axlims = [u.xylims(3) u.xylims(4) u.zmin u.zmax];
    end
    
    if any(isnan(yind))
        yind = 1:m.ny;
    end
    
    zind = 1:nearestpoint(u.zmax*1000,m.cz);
    
    if any(isnan(zind))
        zind = 1:m.nz;
    end
    
    y = [m.y(yind); m.y(yind(end)+1)]/1000;
    z = [m.z(zind); m.z(zind(end)+1)]/1000;
    C = squeeze(log10(m.A(idx,yind,zind)))';
    C = horzcat(C,C(:,end));
    C = vertcat(C,C(end,:));
    
    figure(2);hold on
    [y,z,rho] = plot_cross_section(y,z,C);
    axis(axlims)
    title(['East-West Section @ ',num2str(m.cx(idx)/1000),' km (NS)']);
    set(gca,'Layer','top')
    
    if strcmp(u.gridlines,'off')
        shading flat
    end
    
    %If a data structure exists, project the nearby sites onto the section
    if exist('d','var') && isfield(d,'idy')
        plot(d.y(site_distx<=u.tol)/1000,m.Z(idx,d.idy(site_distx<=u.tol))/1000 -u.zoff,'vk','MarkerFaceColor','k')
    end
    print_figure(['vertical_profiles_',m.niter],['E-W_at_slice_',num2str(idx,'%03.0f')]);

    %Save figure
    figure(1)
    print_figure(['vertical_profiles_',m.niter],['E-W_at_slice_',num2str(idx,'%03.0f'),'_MAP']);
    
    %Output the x,y,z,rho coordinates of the chosen slice 
    x = ones(length(y),1).*m.cx(idx)/1000;
    rho = squeeze(rho);

end

end %END MAIN
