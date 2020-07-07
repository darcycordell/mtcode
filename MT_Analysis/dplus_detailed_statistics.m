function dobs = dplus_detailed_statistics(d)
%Function which calculates the D+ fit for each site in your dataset and
%compiles this into a "dpred" data structure. Also allows you to plot
%misfit statistics and calculate r.m.s. values for D+ fits for each period
%or in map view etc.
%
% Usage: dpred = dplus_detailed_statistics(d)
%
% Inputs: d is a standard data structure
%
% Outputs: dpred is a standard data structure containing the D+ predicted
% data for the impedances, apparent resistivities and phases for the
% off-diagonal components only. (Note: D+ theory does not work for diagonal
% components)
%

[d] = set_map_projection(d);
[L] = load_geoboundary_file_list;

dpred = d;
dpred.Z = nan(size(d.Z));
dpred.Zerr = nan(size(d.Zerr));
dpred.rho = nan(size(d.rho));
dpred.rhoerr = nan(size(d.rhoerr));
dpred.pha = nan(size(d.pha));
dpred.phaerr = nan(size(d.phaerr));
dpred.tip = nan(size(d.tip));
dpred.tiperr = nan(size(d.tiperr));
dpred.name = ['dplus_for_',d.name];

disp('Computing D+ fits for each MT site ...')
tic; obs = []; pred = []; err = [];
for is = 1:d.ns
   
    [detZ, detZerr] = calc_determinant(d.Z(:,:,is),d.Zerr(:,:,is));
    [dpd] = dplus(detZ,detZerr,d.T);
    obs = [obs; real(dpd.Zorig); imag(dpd.Zorig)];
    pred = [pred; real(dpd.Z); imag(dpd.Z)];
    err = [err; real(dpd.Zerr); imag(dpd.Zerr)];
    
    
    [dplus_vars(is,2)] = dplus(d.Z(:,2,is),d.Zerr(:,2,is),d.T);
    [dplus_vars(is,3)] = dplus(-d.Z(:,3,is),d.Zerr(:,3,is),d.T);

    [~,ind] = ismember(dplus_vars(is,2).T,d.T);

    dpred.Z(ind,2,is) = dplus_vars(is,2).Z;
    dpred.rho(ind,2,is) = dplus_vars(is,2).rho;
    dpred.pha(ind,2,is) = dplus_vars(is,2).pha;
    dpred.Zerr(ind,2,is) = dplus_vars(is,2).Zerr;
    
    [~,ind] = ismember(dplus_vars(is,3).T,d.T);
    
    dpred.Z(ind,3,is) = -dplus_vars(is,3).Z;
    dpred.rho(ind,3,is) = dplus_vars(is,3).rho;
    dpred.pha(ind,3,is) = dplus_vars(is,3).pha-180;
    dpred.Zerr(ind,3,is) = dplus_vars(is,3).Zerr;


end
toc
%%
[stats_all] = detailed_statistics(dpred,d);

[stats_det] = statistics(obs,pred,err);

disp(['Overall R.M.S. Misfit = ',num2str(stats_all.rms)])
disp(['Determinant R.M.S. Misfit = ',num2str(stats_det.rms)])
disp('Change D+ error floor in user_defaults to get overall r.m.s. misfit closer to 1.0')

%%
dobs = dpred; %Necessary to flip so that the observed data has error values
dpred = d;


while 1
    
    misfit_menu = menu('','Residual Histograms','Misfit By Period','Misfit Map','Misfit Map (By Period)','Quit');

    if misfit_menu == 1 %Residual histograms-------------------
        close all
        plot_misfit_residual_histograms(dobs,dpred)
        print_figure(['misfit_stats_',dpred.niter],'Residuals'); %Save figure

    elseif misfit_menu == 2 %Misfit by period------------------
        close all
        plot_misfit_rms_by_period(dobs,dpred)

    elseif misfit_menu == 3 %Misfit map view-------------------
        close all;
        plot_misfit_rms_map(dobs,dpred)    
        plot_geoboundaries(L);     
        print_figure(['misfit_stats_',dpred.niter],'rms_map_view'); %Save figure

    elseif misfit_menu == 4 %Misfit map view (by period)

        close all
        plot_misfit_rms_by_period_map(dobs,dpred)

    end

    if misfit_menu == 5
        break
    end
end

end