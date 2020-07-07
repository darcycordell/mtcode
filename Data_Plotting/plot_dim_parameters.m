function plot_dim_parameters(d)
%%
% Function which plots dimensionless parameters such as Swift skew, Bahr
% skew, Swift angle, s1, s2, etc. as a function of period
%
% Usage: plot_dim_parameters(d)
%
% "d" is data structure input in standard format
%
% There is an option to plot strike angle of maximum phase split. However,
% this option takes quite long (5 - 10 s for each station). If it is not
% selected, the plot is left blank.

rot_ang = 0;

u = user_defaults;

dr = rotate_d(d,rot_ang);

phase_split_menu = menu('Calculate Phase Splits? (Takes Longer)','Yes','No');

is = 1;
while 1

    %Calculate dimensionless parameters in "dim" structure
    [dim] = calc_dim_parameters(dr.Z(:,:,is),rot_ang);
    
    if phase_split_menu == 1
        %Calculate max phase split strike angle
        [strike_max_split] = calc_pha_splits(d,is);
    else %Else make NaN
        strike_max_split = nan(d.nf,1);
    end

    %Plot s1, s2, s3, s4
    set_figure_size(1);
    subplot(2,2,1)
    loglog(d.T,abs(dim.s1),'-r'); hold on
    loglog(d.T,abs(dim.s2),'r:') ;
    loglog(d.T,abs(dim.d1),'k:') ;
    loglog(d.T,abs(dim.d2),'k-') ;
    legend('s1','s2','d1','d2')
    ylabel('Dim Parameter Magnitude');
    xlabel('Period (s)')
    title('Dimensionless Parameters')
    axis([u.Tlim(1) u.Tlim(2) 0.5*min(abs([dim.s1; dim.s2; dim.d1; dim.d2])) 2*max(abs([dim.s1; dim.s2; dim.d1; dim.d2]))])

    %Plot Swift Skew and Bahr Skew
    subplot(2,2,2)
    semilogx(d.T,dim.swift_skew,'-r'); hold on
    semilogx(d.T,dim.bahr_skew,'r:') ;
    legend('Swift Skew','Bahr Skew')
    ylabel('Skew');
    xlabel('Period (s)')
    title('Swift and Bahr Skew')
    axis([u.Tlim(1) u.Tlim(2) 0.9*min([dim.swift_skew; dim.bahr_skew]) 1.1*max([dim.swift_skew; dim.bahr_skew])])

    %Plot Swift angle
    subplot(2,2,3)
    semilogx(d.T,dim.alpha,'-r'); hold on
    ylabel('Angle (degrees)');
    xlabel('Period (s)')
    title('Swift Angle')
    axis([u.Tlim(1) u.Tlim(2) 0.9*min(dim.alpha) 1.1*max(dim.alpha)])

    %Plot Strike Angle with maximum phase split (if desired)
    subplot(2,2,4)
    semilogx(d.T,strike_max_split,'-b'); hold on
    ylabel('Strike Angle (deg)');
    xlabel('Period (s)')
    title('Strike Angle with Maximum Phase Split');
    if phase_split_menu == 1
        axis([u.Tlim(1) u.Tlim(2) 0.9*min(strike_max_split) 1.1*max(strike_max_split)])
    end


    annotation('textbox', [0 0.9 1 0.08], ...
        'String', ['Dimensionless Parameters for Site ',d.site{is}], ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center','FontSize',12,'FontWeight','bold')
    
    print_figure(['dim_params_',d.niter],['dim_param_',d.site{is}]); %Save figure

    next_menu = menu('','Next Station','Previous Station','Return');

    if next_menu == 1
        is = is+1;
        if is>d.ns
            is = d.ns;
        end
    elseif next_menu == 2
        is = is-1;
        if is<1
            is = 1;
        end
    else
        break
    end


end
