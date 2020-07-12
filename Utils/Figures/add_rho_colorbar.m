function add_rho_colorbar(u)

    hcb = colorbar;
    set(hcb, 'Ticks', get(hcb,'Ticks')); % need this to ensure number of ticks doesn't change when resizing figure
    if strcmp(u.clabel,'log')    
        hcb.Label.String = 'log10 ( Resistivity \Omega m )';
    elseif strcmp(u.clabel,'linear')
        set(hcb,'TickLabels',num2cell(round(10.^(hcb.Ticks),2,'significant'))')
        hcb.Label.String = 'Resistivity (\Omega m)';
    end

end