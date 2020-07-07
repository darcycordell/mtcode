function make_M3_param_Callback(hObject, ~, ~) % format has changed to allow loading
    H=guidata(hObject);
    
    if isfield(H,'dat')
        filename=['par_',strtok(H.dat,'.'),'.m'];
    else
        warndlg('Need to load data before saving a parameter file.')
        return
    end

    if exist(filename,'file')==2
        filename=uiputfile('par_*.m','Parameter file already exists!',filename);
    end

    fid=fopen(filename,'w');

    fprintf(fid,'%s\n',[get(H.xspacing,'string'),'    % spacing in x (m)']);
    fprintf(fid,'%s\n',[get(H.yspacing,'string'),'     % spacing in y (m)']);
    fprintf(fid,'%s\n',[get(H.moutcX,'string'),'       % mesh out the core in X']);
    fprintf(fid,'%s\n',[get(H.Xinc,'string'),'  % increase in X by']);
    fprintf(fid,'%s\n',[get(H.moutcY,'string'),'       % mesh out the core in Y']);
    fprintf(fid,'%s\n',[get(H.Yinc,'string'),'  % increase in Y by']);
    fprintf(fid,'%s\n',[get(H.firstthick,'string'),' % First Thickness']);
    fprintf(fid,'%s\n',[get(H.nl,'string'),'          % number of layers']);
    fprintf(fid,'%s\n',[get(H.Zinc,'string'),'  % increase thickness by']);
    fprintf(fid,'%s\n',[get(H.airspacing,'string'),'  % air cell thickness']);
    fprintf(fid,'%s\n',[get(H.res,'string'),'  % halfspace resistivity']);
    
    fprintf(fid,'%s\n',[get(H.min_per,'string'),'      % minimum period to use']);
    fprintf(fid,'%s\n',[get(H.max_per,'string'),'  % maximum period to use']);
    fprintf(fid,'%s\n',[get(H.per_skip,'string'),'     % export every nth period']);
    fprintf(fid,'%s\n',[get(H.errflr_Zdiag,'string'),'     % xx and yy error floor (%)']);
    fprintf(fid,'%s\n',[get(H.errflr_Zodiag,'string'),'     % xy and yx error floor (%)']);
    fprintf(fid,'%s\n',[get(H.errflr_tip,'string'),'     %tipper error floor (absolute)']);
    fprintf(fid,'%s\n',[get(H.errflr_type,'string'),'     %error floor type']);
    fprintf(fid,'%s\n',[num2str(get(H.resp,'value')),'     % [1= full tensor] [2= full tensor and tip] [3= off-diag only] [4= tip only]']);
    fprintf(fid,'%s\n',[num2str(get(H.conj_Z,'value')),'     % 0= no conj Z, 1=conj Z']);  
    fprintf(fid,'%s\n',[num2str(get(H.conj_tip,'value')),'     % 0= no conj tip, 1=conj tip']);  

    fclose(fid);
    
    disp([filename,' created'])

    guidata(hObject, H);
    
end