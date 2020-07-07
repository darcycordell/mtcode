function [is] = plot_slice_menus(m,is,d)
%Function to plot model slice-by-slice by going to the next or previous
%slice
%
% Usage: [is] = plot_slice_menus(m,is,d) OR plot_slice_menus(m,is)
%
% "m" is the model structure
% "is" is the index to plot initially
% "d" is an OPTIONAL data structure

if ~exist('d','var')
    d = make_nan_data;
end

irun=1;
while irun==1; %Secondary while loop: Exit loop by selecting slice
clf
plot_slice(m,is,d)

%Option to specify slice by going deeper or shallower. Choose
%the slice which best represents the feature to be edited.
imenu=menu('','Next Slice (Deeper)','Previous Slice (Shallower)','Done');

if imenu==1
        is=is+1;
        if is>=m.nz
            is=m.nz;
            disp('Reached base of model. No deeper slice')
        end
    elseif imenu==2
        is=is-1;
        if is<1
            is=1;
            disp('Reached top of model. No shallower slice')
        end
    else
        return
    end
end

end

