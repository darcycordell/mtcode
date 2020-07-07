function [log] = load_logfile_wsinv(logfile)
%
% Function which loads a WSINV3D *.log file output from the inversion
%
% Usage: [log] = load_logfile_wsinv(logfile)
%
% "log" is a structure containing data fit and model norm info for each iteration
% "logfile" is the ModEM log file name string
%

fid=fopen(logfile,'r');

disp('Loading WSINV Log file')

i=1;
while 1
    
line=fgetl(fid);
if length(line)>11
    if strcmp(line(1:11),'<<< ITER NO')
        line(strfind(line,'=')) = [];
        key = 'RMS';
        idx = strfind(line,key);
        log.rms(i) = sscanf(line(idx(1)+length(key):end),'%g',1);
        
        key = 'LM';
        idx = strfind(line,key);
        log.m2(i) = sscanf(line(idx(1)+length(key):end),'%g',1);
        
        i=i+1;


    end
end

if line == -1
    break
end

end

log.niter = length(log.rms);

%These variables do not exist in the WSINV logfile and are set to NaN. They
%only exist in the ModEM inversion logfile.
log.f = nan(log.niter,1);
log.lambda = nan(log.niter,1);
log.alpha = nan(log.niter,1);

end
%====================================================================