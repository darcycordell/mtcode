function [log] = load_logfile_modem(logfile)
%
% Function which loads a ModEM *.log file output from the inversion
%
% Usage: [log] = load_logfile_modem(logfile)
%
% "log" is a structure containing data fit and model norm info for each iteration
% "logfile" is the ModEM log file name string
%

fid=fopen(logfile,'r');

disp('Loading ModEM Log file')

line=fgetl(fid);
line_num=1;
initial_flag=true;
    while line ~= -1

        % first read in the initial inversion parameters
        if initial_flag
            ind_tmp=strfind(line,'START:');
            if ~isempty(ind_tmp)
                ind=findstr(line,'f=');
                if ~isempty(ind)
                    f_init=sscanf(line((ind+2):(ind+14)),'%E');
                else
                    f_init=NaN;
                end

                ind=findstr(line,'m2=');
                if ~isempty(ind)
                    m2_init=sscanf(line((ind+3):(ind+15)),'%E');
                else
                    m2_init=NaN;
                end

                ind=findstr(line,'rms=');
                if ~isempty(ind)
                    rms_init=sscanf(line((ind+4):(ind+16)),'%f');
                else
                    rms_init=NaN;
                end

                ind=findstr(line,'lambda=');
                if ~isempty(ind)
                    lambda_init=sscanf(line((ind+7):(ind+19)),'%E');
                else
                    lambda_init=NaN;
                end

                ind=findstr(line,'alpha=');
                if ~isempty(ind)
                    alpha_init=sscanf(line((ind+6):(ind+17)),'%E');
                else
                    alpha_init=NaN;
                end

                initial_flag=false; % make sure the initial conditions do not get read in again
            end
        end

        % now read in the results of each iteration
        ind_tmp=strfind(line,'Completed NLCG iteration'); %find the line where the results from each iteration are written
        if ~isempty(ind_tmp)
            iter_num=sscanf(line(ind_tmp+24:length(line)),'%i');
            line=fgetl(fid);
            % if more than one node is used in the inversion (ran as parallel),
            % it seems like the output lines in the logfile get mixed up and this following
            % string sometimes gets written in over the data we are trying to
            % read, making it fail.
            ind_decrease=strfind(line,'Sufficient decrease condition satisfied, exiting line search');
            if isempty(ind_decrease)

                ind=findstr(line,'f=');
                if ~isempty(ind)
                    f(iter_num)=sscanf(line((ind+2):(ind+14)),'%E');
                else
                    f(iter_num)=NaN;
                end

                ind=findstr(line,'m2=');
                if ~isempty(ind)
                    m2(iter_num)=sscanf(line((ind+3):(ind+15)),'%E');
                else
                    m2(iter_num)=NaN;
                end

                ind=findstr(line,'rms=');
                if ~isempty(ind)
                    rms(iter_num)=sscanf(line((ind+4):(ind+16)),'%f');
                else
                    rms(iter_num)=NaN;
                end

                ind=findstr(line,'lambda=');
                if ~isempty(ind)
                    lambda(iter_num)=sscanf(line((ind+7):(ind+19)),'%E');
                else
                    lambda(iter_num)=NaN;
                end

                ind=findstr(line,'alpha=');
                if ~isempty(ind)
                    if line(ind+6)=='*';% sometimes the alpha value is recorded as '************', we will skip these values
                        alpha(iter_num)=NaN;
                    else
                        alpha(iter_num)=sscanf(line((ind+6):(ind+17)),'%E');
                    end
                else
                    alpha(iter_num)=NaN;
                end

            else %if the 'sufficient decrease....' line overwrote the data we are trying to read in, just make them NaN
                f(iter_num)=NaN;
                m2(iter_num)=NaN;
                rms(iter_num)=NaN;
                lambda(iter_num)=NaN;
                alpha(iter_num)=NaN;
            end

        end

        line=fgetl(fid); %read in the next line
        line_num=line_num+1; %update the line number (not used for anything but debugging)



    end
    
    if exist('f','var') == 1
        log.f = f';
        log.m2 = m2';
        log.rms = rms';
        log.lambda = lambda';
        log.alpha = alpha';
        log.niter = length(f);
    else
        log.f = f_init';
        log.m2 = m2_init';
        log.rms = rms_init';
        log.lambda = lambda_init';
        log.alpha = alpha_init';
        log.niter = length(f_init);
    end
end
%====================================================================