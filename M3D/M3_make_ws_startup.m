function M3_make_ws_startup(hObject, ~, ~)
% Writes startup file for wsinv3dmt.    

    H=guidata(hObject);
          
    [data_file,~] = uigetfile({'*'},' choose data file');
    [model_file,~] = uigetfile({'*'},' choose initial model file');

    s_prompt = {'data file:','initial model file:','output file:','priori model file(leave empty for default):','control model index(leave empty for default):','Target R.M.S.:','Max number of iterations:'...
        ,'Model length scale(leave empty for default):','Lagrange info(leave empty for default):','Error tol. level(leave empty for default):', 'Number of responses (12, 8, or 4)','Z error floor','Tipper error floor'};
    s_titles  = 'Choose startup parameters';
    out_file=['out_',model_file];def='default';
    s_def= {data_file,model_file,out_file,'','','1.','10','','','','12','5','15'}; %Turkiye
    def_lim = char(inputdlg(s_prompt,s_titles,1,s_def));
    for spa=4:10
        if def_lim(spa,1:2)=='  '
            def_lim(spa,1:7)='default';
        end
    end
    nr=str2double(def_lim(11,:));
    if nr == 12
        inv_type='5';
    elseif nr== 8
        inv_type='1';
    elseif nr == 4
        inv_type='2';
    end %there are also only tipper; and Zxy, Zyx, tipper inversions possible as well (define yourself!)
    line{1}=['INVERSION_TYPE      ',inv_type];
    line{2}=['DATA_FILE           ',def_lim(1,:)];
    line{3}=['MIN_ERROR_Z         ',def_lim(12,:)];
    line{4}=['MIN_ERROR_T         ',def_lim(13,:)];
    line{5}=['OUTPUT_FILE         ',def_lim(3,:)];
    line{6}=['INITIAL_MODEL_FILE  ',def_lim(2,:)];
    line{7}=['PRIOR_MODEL_FILE    ',def_lim(4,:)];
    line{8}=['CONTROL_MODEL_INDEX ',def_lim(5,:)];
    line{9}=['TARGET_RMS          ',def_lim(6,:)];
    line{10}=['MAX_NO_ITERATION    ',def_lim(7,:)];
    line{11}=['MODEL_LENGTH_SCALE  ',def_lim(8,:)];
    line{12}=['LAGRANGE_INFO       ',def_lim(9,:)];
    line{13}=['ERROR_TOL_LEVEL     ',def_lim(10,:)];
    fid=fopen('startup','w');
    if nr ==12
        for poi=1:13
            fprintf(fid,'%s\n',line{poi});
        end
    else %no tipper error is output
        for poi=1:3
            fprintf(fid,'%s\n',line{poi});
        end
        for poi=4:13
            fprintf(fid,'%s\n',line{poi});
        end
    end
    
    fclose(fid);      
    guidata(hObject, H);
    
end