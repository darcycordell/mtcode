function plot_mesh(m,d)
% Function to plot the mesh without using pcolor. This is useful to check that
% stations are located in the cell center
%
% Usage: plot_mesh(m,d)
%
%

set_figure_size(1);
plot_mesh_option(m,d,1);
plot_mesh_option(m,d,0);
manual_legend('Cell Centers','.k','Cell Edges','-g','Station Locations','vr');

end

function plot_mesh_option(m,d,flag)

u = user_defaults;
[L] = load_geoboundary_file_list;

if flag
    subplot(1,2,1);
else
    subplot(1,2,2);
end

plot(m.Xc(:)/1000,m.Yc(:)/1000,'.k'); hold on

for i = 1:m.nx+1
    plot(m.X(i,:)/1000,m.Y(i,:)/1000,'-g');
end

for i = 1:m.ny+1
    plot(m.X(:,i)/1000,m.Y(:,i)/1000,'-g');
end

plot(d.x/1000,d.y/1000,'vr','MarkerFaceColor','r')
plot_geoboundaries(L,d.origin,0)

axis equal

if flag
    %Set plotting limits
    if length(u.xylims)==4 %If xy limits were specified in user_defaults
        xind = nearestpoint(u.xylims(1),m.cx/1000,'previous'):nearestpoint(u.xylims(2),m.cx/1000,'next');
        yind = nearestpoint(u.xylims(3),m.cy/1000,'previous'):nearestpoint(u.xylims(4),m.cy/1000,'next');
        axlims = [u.xylims(3) u.xylims(4) u.xylims(1) u.xylims(2)];

    else %If the limits specified are not a vector of numbers with length 4 then just take non-padding cells
        xind = m.npad(1)+1:m.nx-m.npad(1);
        yind = m.npad(2)+1:m.ny-m.npad(2);
        axlims = [sort([m.y(min(yind)) m.y(max(yind))]) sort([m.x(min(xind)) m.x(max(xind))])]/1000;
    end

    axis(axlims)
    title('Map of Mesh Survey Area');
else
    axis([m.y(1) m.y(end) m.x(1) m.x(end)]/1000)
    title('Map of Full Mesh Area')
end
set(gca,'Layer','top')
xlabel('Distance East-West (km)')
ylabel('Distance North-South (km)')

end

