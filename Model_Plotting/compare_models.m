function [a,b,da,db] = compare_models(a,b,da,db,flag)
% Function which allows for the comparison of two models which are both in
% a standard model structure. Interpolation is carried out over the area
% shared by both meshes. Plots of each model are made along with a 2D
% histogram to compare models. Pearson correlation coefficient is
% calculated to compare the similarity between models. 
%
% Usage: [a,b,da,db] = compare_models(a,b,da,db)
%
% "a" is the first model structure (Model #1)
% "b" is the second model structure (Model #2)
% "da" is the data structure associated with the first model structure
% "db" is the data structure associated with the second model structure
% "flag" is to plot (1) or not plot (0). Plotting is default.
%
% The outputs (a,b,da,db) are the the interpolated versions of the inputs
%
% The data structures (da,db) are necessary to geo-referenced the center of
% the meshes in latitude and longitude space
%
% Note: Both meshes MUST be referenced to the same point in z-space. For
% example, both must be referenced to sea level.
%
% For debugging purposes it is often useful to keep the original model and
% data structures before interpolation
% a_orig = a; da_orig = da;
% b_orig = b; db_orig = db;
%
%
% (C) 2020 Unsworth Research Group (University of Alberta, Edmonton, Canada)

if ~exist('flag','var')
    flag = 1;
end

%Determine the NS and EW distance offset in meters between the two mesh origins
dist_ns=distance('gc',da.origin(1),da.origin(2),db.origin(1),da.origin(2),6371000);
dist_ew = distance('gc',da.origin(1),da.origin(2),da.origin(1),db.origin(2),6371000);

%Move the b mesh origin to match the a mesh origin.
% In order to do this, we need to determine if the b mesh origin is NW, NE,
% SW or SE of the a mesh origin
if db.origin(1) >= da.origin(1) && db.origin(2) >= da.origin(2)
    %b mesh origin is NE of "a" mesh origin
    ew = -dist_ew;
    ns = -dist_ns;
elseif db.origin(1) > da.origin(1) && db.origin(2) < da.origin(2)
    %b mesh origin is NW of "a" mesh origin
    ew = dist_ew;
    ns = -dist_ns;
elseif db.origin(1) <=da.origin(1) && db.origin(2) <= da.origin(2)
    %b mesh origin is SE of "a" mesh origin
    ew = dist_ew;
    ns = dist_ns;
elseif db.origin(1)<=da.origin(1) && db.origin(2) > da.origin(2)
    %b mesh origin is SW of "a" mesh origin
    ew = -dist_ew;
    ns = dist_ns;    
end

% Move the b mesh vectors to align with the "a" mesh and set origins
% accordingly
b.cx = b.cx-ns;
b.cy = b.cy-ew;
db.origin = da.origin;
db.x = db.x-ns;
db.y = db.y-ew;

% Currently, models cannot be compared if one has topography and the other
% does not
if b.origin(3) == 0 && a.origin(3) ~=0 || b.origin(3) ~=0 && a.origin(3) == 0
    disp('One of your models has topography while the other does not. Make sure your models are geo-referenced properly!')           
end
  
interp_menu = menu('','Interpolate Onto 1st Model Mesh','Interpolate Onto 2nd Model Mesh');

% Find the shared region that the two models share
if interp_menu == 1 %a.cx(end)-a.cx(1) < b.cx(end)-b.cx(1)
    %Use "a" mesh NS range
    cx = a.cx; b.npad(1) = a.npad(1); x = a.x;
    indx = a.npad(2)+1:a.nx-a.npad(2); %x range to plot model (ignore padding)
else
    %Use "b" mesh NS range
    cx = b.cx; a.npad(1) = b.npad(1); x = b.x;
    indx = b.npad(2)+1:b.nx-b.npad(2);
end

if interp_menu == 1 %a.cy(end)-a.cy(1) < b.cy(end)-b.cy(1)
    %Use "a" mesh EW range
    cy = a.cy; b.npad(2) = a.npad(2); y = a.y;
    indy = a.npad(2)+1:a.ny-a.npad(2); %y range to plot model (ignore padding)
else
    %Use "b" mesh EW range
    cy = b.cy; a.npad(2) = b.npad(2); y = b.y;
    indy = b.npad(2)+1:b.ny-b.npad(2);
end

if interp_menu == 1 %abs(a.cz(end)-a.cz(1)) < abs(b.cz(end)-b.cz(1))
    %Use "a" mesh vertical range
    cz = a.cz; z = a.z;
    dx = a.dx; dy = a.dy; dz = a.dz;
    Z = a.Z;
else
    %Use "b" mesh vertical range
    cz = b.cz; z = b.z;
    dx = b.dx; dy = b.dy; dz = b.dz;
    Z = b.Z;
end

if interp_menu == 1
    origin = a.origin;
else
    origin = b.origin;
end

nx = length(cx);
ny = length(cy);
nz = length(cz);

%Create a 3D meshgrid of shared model space
[Xq, Yq, Zq] = meshgrid(cy,cx,cz);

%Interpolate each mesh onto the new model space
Aq = interp3(a.cy,a.cx,a.cz,a.A,Xq,Yq,Zq);
Bq = interp3(b.cy, b.cx, b.cz, b.A, Xq, Yq, Zq);

%Reset model structures with new variables
a.A = Aq;
a.cx = cx; a.cy = cy; a.cz = cz;
a.x = x; a.y = y; a.z = z;
a.nx = nx; a.ny = ny; a.nz = nz;
a.dx = dx; a.dy = dy; a.dz = dz;
a.origin = origin;
a.Z = Z;

b.niter = '';

b.A = Bq;
b.cx = cx; b.cy = cy; b.cz = cz;
b.x = x; b.y = y; b.z = z;
b.nx = nx; b.ny = ny; b.nz = nz;
b.dx = dx; b.dy = dy; b.dz = dz;
b.origin = origin;
b.Z = Z;

[a.X,a.Y] = meshgrid(cy,cx);
b.X = a.X;
b.Y = a.Y;

%Done Interpolation--------------------------------------------------------


%Begin plotting sequence---------------------------------------------------
if flag
iz = round(b.nz/2); %z slice to plot

while 1
    
A = a.A(indx,indy,iz);
B = b.A(indx,indy,iz);

if size(A) == size(B) %If A ~= B, something has gone wrong with interpolation!

    set_figure_size(1);
    
    %Plot Model #1
    ax2 = subplot(1,3,1);
    cla(ax2);
    plot_slice(a,iz,da);
    title(['Model #1: Depth=',num2str(cz(iz)./1000),' km. Slice #',num2str(iz) ]);
    
    %Plot Model #2
    ax3 = subplot(1,3,2);
    cla(ax3);
    plot_slice(b,iz,db);
    
    title(['Model #2: Depth=',num2str(cz(iz)./1000),' km. Slice #',num2str(iz)]);

    
    %Linearize models into a single vector
    Avec = A(:);
    Bvec = B(:);

    %Create a 2D histogram
    vXEdge = linspace(-2,4,50);
    vYEdge = linspace(-2,4,50);
    
    %FOR SEISMIC
    %vXEdge = linspace(2,4,50);
    
    mHist2d = hist2d(log10(Avec),log10(Bvec),vYEdge,vXEdge);  

    nXBins = length(vXEdge);
    nYBins = length(vYEdge);
    vXLabel = 0.5*(vXEdge(1:(nXBins-1))+vXEdge(2:nXBins));
    vYLabel = 0.5*(vYEdge(1:(nYBins-1))+vYEdge(2:nYBins));
    
    %Plot 2D Histogram
    ax1 = subplot(1,3,3);
    colormap(ax1,flip(gray)); colorbar;     
    pcolor((vXLabel), (vYLabel),log10(mHist2d)); shading flat
    xlabel('Model #2');ylabel('Model #1'); 
    set(gca,'dataaspectratio',[1,1,1])
    
    if all(isnan(Avec))==1 || all(isnan(Bvec)) == 1
        %If you are in only topography cells and either of the models are
        %all NaN do not do correlation
        
    else %Calculate correlation coefficents
        
        %Lots of options for doing correlation. Pearson seems to be best
        
        d1 = (Avec);
        d2 = (Bvec);

        d1(isnan(d2)) = [];
        d2(isnan(d2)) = [];   

        d2(isnan(d1)) = [];
        d1(isnan(d1)) = [];

        %Correlation Coefficient and R-square are BAD methods of measuring the
        %sameness of models. If both models are ~99 Ohm m, it will give very
        %low correlations!
        %p1 = polyfit(d1,d2,1); % y = M, x = Z
        %fit1 = polyval(p1,d1);
        %resid1 = d2 - fit1;
        %SSresid1 = sum(resid1.^2);
        %SStotal1 = (length(d2)-1)*var(d2);
        %rsq1 = 1-(SSresid1/SStotal1);

        %p2 = polyfit(d2,d1,1); % y = Z, x = M
        %fit2 = polyval(p2,d2);
        %resid2 = d1 - fit2;
        %SSresid2 = sum(resid2.^2);
        %SStotal2 = (length(d1)-1)*var(d1);
        %rsq2 = 1-(SSresid2/SStotal2);

        %Try Cosine Similarity instead (0 - 1)
        %cos_sim = sum(d1.*d2)/(sqrt(sum(d1.^2))*sqrt(sum(d2.^2)));

        %Try Pearson Correlation Coefficient (-1 to +1)
        pearson = sum((d1-mean(d1)).*(d2-mean(d2)))/(sqrt(sum((d1-mean(d1)).^2)*sum((d2-mean(d2)).^2)));

        %title(['Cosine Similarity = ',num2str(cos_sim)])%----------------<
        title(['Pearson Correlation = ',num2str(pearson)])
    end
  
    %Uncomment to print figure
    print_figure(['compare_model_',da.niter],['compare_model_',num2str(iz)])
    
    main_menu = menu('','Next Slice','Previous Slice','Choose Area to Correlate','Quit');

    if main_menu == 1 %View next slice
        iz = iz+1;
        if iz > nz
            iz = nz;
        end
        clf

    elseif main_menu == 2 %View previous slice
        iz = iz-1;
        if iz < 1
            iz = 1;
        end
        clf

    elseif main_menu ==3 %Click outline to correlate
        %You can click on either of the models to find the correlation
        %between a particular area of interest.
        
        disp('Click outline of area to paint on horizontal slice. Right click to stop selection.')

        figure(1);
        [hx,hy]=ginputExtra(100); %select polygons on horizontal slice

        disp(['Points clicked .... x points: [', num2str(hx),'], y points: [',num2str(hy),']']);

        if isempty(hx) == 1 || length(hx) < 3
            error('Polygon must be defined by at least 3 points!')
        end

        %Find the indices of the A matrix which are contained within the
        %polygon
        [X,Y]=meshgrid(a.y,a.x);

        [in] = inpolygon(X/1000,Y/1000,hx,hy);
        [row, col]=find(in==1);

        %Set new x and y indices to correlate
        indx = min(row):max(row);
        indy = min(col):max(col); 

    else
        %close all
        break
    end

end

end %END While Loop

end 

end %END MAIN