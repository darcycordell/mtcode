function [invtype, dobs_file, Z_err_floor, tip_err_floor, logfile, m0_file] = load_startup_wsinv(startupfile)

fid = fopen(startupfile,'r');
tmp=0;
Z_err_floor = 0;
tip_err_floor = 0;

while isempty(tmp)~=1
    tmp = fscanf(fid,'%s',1); % reading
    if strcmp(tmp,'INVERSION_TYPE')
        invtype = str2double(fscanf(fid,'%s',1));
    elseif strcmp(tmp,'DATA_FILE')
        dobs_file = fscanf(fid,'%s',1);
    elseif strcmp(tmp,'MIN_ERROR_Z')
        Z_err_floor = str2double(fscanf(fid,'%s',1));
    elseif strcmp(tmp,'MIN_ERROR_T')
        tip_err_floor = str2double(fscanf(fid,'%s',1));
        
    elseif strcmp(tmp,'OUTPUT_FILE')
        logfile = [fscanf(fid,'%s',1),'.log'];
    elseif strcmp(tmp,'INITIAL_MODEL_FILE')
        m0_file = fscanf(fid,'%s',1);
    end
end


end