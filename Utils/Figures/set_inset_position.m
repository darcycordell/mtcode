function set_inset_position(ax,h_inset,loc)
% Set Rose diagram inset location (top left, top right, bottom left or 
% bottom right)
% Inputs:
% ax = position of the main figure, retrieved before function call with ax=get(main_fig,'Position');
% h_inset = handle of the axes for the inset figure
% loc = string or xy coordinate location of inset (from user_defaults)
%
% this function is currently only called by plot_phase_tensor_map and
% plot_induction_vector_map
        
u = user_defaults;
inset_size=u.inset_size*.7;

if ~isnumeric(loc)
    switch loc
        case 'SW'
            set(h_inset,'Position', [ax(1) ax(2) inset_size inset_size]) % SW    
        case 'SE'
            set(h_inset,'Position', [ax(1)+ax(3)-inset_size ax(2) inset_size inset_size]) % SE
        case 'NE'
            set(h_inset,'Position', [ax(1)+ax(3)-inset_size ax(2)+ax(4)-inset_size inset_size inset_size]) % NE
        case 'NW'
            set(h_inset,'Position', [ax(1) ax(2)+ax(4)-inset_size inset_size inset_size]) % NW
        otherwise
            disp('Unrecognized format of u.inset_loc. Check your user_defaults file')
            set(h_inset,'Position', [ax(1) ax(2)+ax(4)-inset_size inset_size inset_size]) % NW
    end
elseif isnumeric(loc) && numel(loc)==2
    set(h_inset,'Position', [loc(1) loc(2) inset_size inset_size]) % custom
else
    disp('Unrecognized format of u.inset_loc. Check your user_defaults file')
    set(h_inset,'Position', [ax(1) ax(2)+ax(4)-inset_size inset_size inset_size]) % NW
end
        
end