function [log] = load_logfile_2DNLCG(logfile)
%
% Function which loads a ModEM *.log file output from the inversion
%
% Usage: [log] = load_logfile_modem(logfile)
%
% "log" is a structure containing data fit and model norm info for each iteration
% "logfile" is the ModEM log file name string
%

fid=fopen(logfile,'r');

disp('Loading 2DNLCG Log file')

line = fgetl(fid);

while ~strcmp(line(2:6),'iter ')
    line = fgetl(fid);
end

i = 1; num = [];
while 1
    line = fgetl(fid);
    if line ~= -1
        num(i,:) = str2num(line);
        if any(isinf(num(i,:)))
            f = strfind(line,'E');
            line = [line(1:f(1)+3),' ',line(f(1)+4:end)];
            num(i,:) = str2num(line);
        end
        i = i+1;
    else
        break
    end
    
end

N = length(num(:,1));
log.f = num(:,7);
log.m2 = nan(N,1);
log.rms = num(:,3);
log.lambda = nan(N,1);
log.alpha = nan(N,1);
log.niter = N;

end
