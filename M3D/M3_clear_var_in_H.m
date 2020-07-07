function clear_variables_in_H(hObject,H)

    if isfield(H,'x')
        H = rmfield ( H,'x');
    end
    if isfield(H,'y')
        H = rmfield ( H,'y');
    end
    if isfield(H,'per')
        H = rmfield ( H,'per');
    end
    if isfield(H,'nx')
        H = rmfield ( H,'nx');
    end
    if isfield(H,'ny')
        H = rmfield ( H,'ny');
    end
    if isfield(H,'nz')
        H = rmfield ( H,'nz');
    end
    if isfield(H,'Z')
        H = rmfield ( H,'Z');
    end
    if isfield(H,'XX')
        H = rmfield ( H,'XX');
    end
    if isfield(H,'YY')
        H = rmfield ( H,'YY');
    end
    if isfield(H,'x_orig')
        H = rmfield ( H,'x_orig');
    end
    if isfield(H,'y_orig')
        H = rmfield ( H,'y_orig');
    end
    if isfield(H,'lay')
        H = rmfield ( H,'lay');
    end
    if isfield(H,'dat')
        H = rmfield ( H,'dat');
    end
    if isfield(H,'AA')
        H = rmfield ( H,'AA');
    end
    if isfield(H,'mesh_rot')
        H = rmfield ( H,'mesh_rot');
    end
    if isfield(H,'res')
        H = rmfield ( H,'res');
    end
    if isfield(H,'use_it')
        H = rmfield ( H,'use_it');
    end
    if isfield(H,'lay_to_c')
        H = rmfield ( H,'lay_to_c');
    end
    if isfield(H,'f_t')
        H = rmfield ( H,'f_t');
    end

    guidata(hObject, H);
end % end clear_variables_in_H