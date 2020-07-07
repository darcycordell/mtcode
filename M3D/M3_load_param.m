function load_M3_param_Callback(hObject,~,~) % works with new par file format
    
    H=guidata(hObject);

    pause(0.1)
    [parname,path]=uigetfile('par_*.m','Read in parameter file');
    par = load([path,parname]);

    set(H.xspacing,'string',num2str(par(1)));
    set(H.yspacing,'string',num2str(par(2)));     % spacing in y (m)
    set(H.moutcX,'string',num2str(par(3)));       % mesh out the core in X
    set(H.Xinc,'string',num2str(par(4)));  % increase in X by
    set(H.moutcY,'string',num2str(par(5)));       % mesh out the core in Y
    set(H.Yinc,'string',num2str(par(6)));  % increase in Y by
    set(H.firstthick,'string',num2str(par(7))); % First Thickness
    set(H.nl,'string',num2str(par(8)));          % number of layers
    set(H.Zinc,'string',num2str(par(9)));  % increase thickness by
    set(H.airspacing,'string',num2str(par(10)));    % air cell thickness
    set(H.res,'string',num2str(par(11)));   % halfspace resistivity
    set(H.min_per,'string',num2str(par(12)));      % minimum period to use
    set(H.max_per,'string',num2str(par(13)));  % maximum period to use
    set(H.per_skip,'string',num2str(par(14)));     % periods to skip
    set(H.errflr_Zdiag,'string',num2str(par(15))); % xx and yy error floor (%)
    set(H.errflr_Zodiag,'string',num2str(par(16))); % xy and yx error floor (%)
    set(H.errflr_tip,'string',num2str(par(17)));    %tipper error floor (absolute)
    set(H.errflr_type,'string',num2str(par(18)));
    set(H.resp,'value',par(19));     % 1=FT only, 2=FT+tip, 3=off diag Z only, 4=tip only
    set(H.conj_Z,'value',par(20));      % 0 = don't conjugate Z, 1= conjugate Z
    set(H.conj_tip,'value',par(21));      % 0 = don't conjugate tipper, 1= conjugate tipper
    
    disp(['File ',parname,' successfully loaded.'])

    guidata(hObject, H);

end % end load_M3_param_Callback