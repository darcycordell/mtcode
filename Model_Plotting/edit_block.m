function [A,R]=edit_block(R,iminx,imaxx,iminy,imaxy,iminz,imaxz,A)
% Function which takes x,y, and z index ranges and replaces that block with
% a specfied resistivity value within a matrix "A". Also the option to only
% replace conductors with a resistivity less than some value or replace
% resistors with a resistivity greater than some value within the block
% (this option can be slow if the block is large).
%
% Usage: [A,R]=edit_block(R,iminx,imaxx,iminy,imaxy,iminz,imaxz,A)
%
% "A" is the model of resistivities in nx x ny x nz
% "R" is the list of indices which have been replaced
% iminx, imaxx, iminy, imaxy, iminz and imaxz are the x, y, and z indices
% which bound the block to be replaced

 %Option to edit conductors, resistors or replace rho directly
iopt = menu('','Edit conductors','Edit resistors','Replace rho directly');

if iopt==1 %----------> Edit conductors
    %Everything within the block with a resistivity value less than
    %the specified resistivity will be set to the new value
    prompt={'Edit out conductors in box with rho <=','Replace with rho ='};
    dlg_title='Replace rho';
    def={'30','30'};
    num_lines=1;
    new_rho = inputdlg(prompt,dlg_title,num_lines,def);
    %Edit box if model is less than given value
    for p=iminx:imaxx
        for q=iminy:imaxy
            for r=iminz:imaxz
                if A(p,q,r)<=str2double(new_rho{1})
                    if A(p,q,r)<10^15
                        A(p,q,r)=str2double(new_rho{2});

                        %Determine the indices which were replaced and add them and their
                        %resistivity values to R
                        ind = sub2ind(size(A),p,q,r);
                        R(ind,:) = [1; str2double(new_rho{2})];
                    end
                end
            end
        end
    end
    fprintf(['If the resistivity value within a cell within the drawn box is less than ',new_rho{1},' Ohm m,\nthose cells will be rewritten to', new_rho{2}, ', Ohm m otherwise the model is left unchanged.']);

elseif iopt==2 %---------> Edit resistors
    %Everything within the block with a resistivity value greater than
    %the specified resistivity will be set to the new value
    prompt={'Edit out resistors in box with rho >=','Replace with rho ='};
    dlg_title='Replace rho';
    def={'30','30'};
    num_lines=1;
    new_rho = inputdlg(prompt,dlg_title,num_lines,def);
    %Edit box if model is greater than given value
    for p=iminx:imaxx
        for q=iminy:imaxy
            for r=iminz:imaxz
                if A(p,q,r)>=str2double(new_rho{1})
                    if A(p,q,r)<10^15
                        A(p,q,r)=str2double(new_rho{2});

                        %Determine the indices which were replaced and add them and their
                        %resistivity values to R
                        ind = sub2ind(size(A),p,q,r);
                        R(ind,:) = [1; str2double(new_rho{2})];
                    end
                end
            end
        end
    end

   fprintf(['If the resistivity value within a cell within the drawn box is greater than ',new_rho{1},' Ohm m,\nthose cells will be rewritten to', new_rho{2}, ', Ohm m otherwise the model is left unchanged.']);

else %----------> Replace all values in the selected volume
    prompt={'Replace rho values in box with rho = '};
    dlg_title='Replace rho';
    def={'30'};
    num_lines=1;
    new_rho = inputdlg(prompt,dlg_title,num_lines,def);
    
    if isempty(new_rho)
        return
    end

    %Replace rho value directly for the range in a box
    A(iminx:imaxx,iminy:imaxy,iminz:imaxz)=str2double(new_rho{1});

    %Determine the indices which were replaced and add them and their
    %resistivity values to R
    r = reshape(R(:,1),size(A));
    r(iminx:imaxx,iminy:imaxy,iminz:imaxz) = 1;
    R(:,1) = r(:);

    r = reshape(R(:,2),size(A));
    r(iminx:imaxx,iminy:imaxy,iminz:imaxz) = str2double(new_rho{1});
    R(:,2) = r(:);

    fprintf(['All cells within the drawn box are replaced with ',new_rho{1},'.']);
    
end
    
end

