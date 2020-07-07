function edit_covariance(m)
% Function which allows you to "edit" the covariance file for a given model
% The function allows you to view the covariance matrix slice-by-slice
%
% Usage: edit_covariance(m)
%
% Inputs: "m" is a standard model structure
%
%
% To "edit" the covariance file you need to have a *.mat file containing
% the indices of the model which you want to edit in the format ix,iy,iz
% ix are the x indices to edit, iy are the y indices, and iz are the z
% indices. The indices are problem-specific so a separate script needs to
% be written to get the indices you want. One option is to use edit_model.m
% to get the indices (they are saved automatically after editing a model)
%
% The covariance file contains an integer value for each model cell:
%
% 0 = air cells (reserved)
% 1 = normal model cells (reserved)
% 2 - 8 = exceptions
% 9 = fixed cells (reserved; e.g. ocean cells)
%
% You can choose to fix cells by setting them to 9
% You can choose up to 7 "exceptions" (2-8) to specify certain cells with
% different areas
%
% After you make the edited covariance you can save it. Then, you must
% MANUALLY edit the covariance file to specify which exceptions you want
%
% For example if you want to turn off smoothing between 1 and 2 then on
% Line #6 of the covariance file you set it to 1 (# of exceptions)
% and below that on Line #7 you say 1 2 0 (to turn off smoothing between 1
% and 2)
%
cov = ones(m.nx,m.ny,m.nz);
cov(isnan(m.A)) = 0; %Set air cells to zero

a = m;
a.A = cov*10;

exception_count = 2;
is = round(m.nz/2);
while 1
set_figure_size(1);
plot_slice(a,is);
colorbar off
caxis([0 2])
h = manual_legend('Halfspace','sg','Air Cells','sr','Exceptions','sb');

set(h(1),'MarkerFaceColor','g')
set(h(2),'MarkerFaceColor','r')
set(h(3),'MarkerFaceColor','b')

imenu=menu('','Next Slice (Deeper)','Previous Slice (Shallower)','Load Exception Indices (*.mat)','Load Fixed Cell Indices (*.mat)','Save Covariance File','Done');

if imenu==1
    is=is+1;
    if is>=a.nz
        is=a.nz;
        disp('Reached base of model. No deeper slice')
    end
elseif imenu==2
    is=is-1;
    if is<1
        is=1;
        disp('Reached top of model. No shallower slice')
    end
    
elseif imenu==3
    
    [idx_file]=uigetfile({'*.mat'},'Load File Containing Model Exception Indices'); if idx_file == 0; return; end
    
    load(idx_file);
    
    indices = sub2ind(size(m.A),ix,iy,iz);
    
    cov(indices) = exception_count; %Set "exception" cells
    
    a.A = cov*10;
    
    exception_count = exception_count+1;
    
elseif imenu==4
  
    [idx_file]=uigetfile({'*.mat'},'Load File Containing Model Fixed Indices'); if idx_file == 0; return; end
    
    load(idx_file);
    indices = sub2ind(size(m.A),ix,iy,iz);
    
    cov(indices) = 9; %Set "exception" cells to 9 (fixed) or 2 (smoothing exception)
    
    a.A = cov*10;
    
elseif imenu==5
    
    prompt = {sprintf('Enter .cov file settings\n\nFile name'),'x smoothing (Default: 0.3)','y smoothing (Default: 0.3)','z smoothing (Default: 0.3)',sprintf('number of times smoothing applied\n(Default: 1)')};
    dlg_title = 'Cov File';
    def = {'modem_edit.cov','0.3','0.3','0.3','1'};
    answer = inputdlg(prompt,dlg_title,1,def);
    
    cfile=answer{1};
    xsmooth = str2double(answer{2})*ones(m.nz,1); % these need to be vectors
    ysmooth = str2double(answer{3})*ones(m.nz,1);
    zsmooth = str2double(answer{4});
    nsmooth = str2double(answer{5});
    
    write_modem_covariance(cfile,cov,xsmooth,ysmooth,zsmooth,nsmooth);
    
else
    return
end


end