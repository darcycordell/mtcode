function write_model_wsinv_test(outputfile,nx,ny,nz,dx,dy,dz,A,niter,RMS)
%Function saves model in wsinv format which can be then read into WSINV
%inversion
%
% Usage: write_model_wsinv(outputfile,nx,ny,nz,dx,dy,dz,A,niter,RMS)
%
% Inputs:
%       outputfile: Output filename string
%       nx, ny, nz are the number of cells in the NS, EW and Z directions,
%       respectively
%
%       dx,dy,dz are the vectors of cell thicknesses (e.g. dx should be a
%       nx x 1 vector of cell thicknesses in meters)
%
%       A is an nx x ny x nz matrix of resistivity values
%
%       niter is the iteration umer
%       RMS is the rms of that model
%
%
    fid=fopen(outputfile,'wt');
    fprintf(fid,'%s','#Iteration No.'); % mjc, no first space
    fprintf(fid,'%4.0f',niter);
    fprintf(fid,'%s','  RMS =       '); % mjc, rms from 1b
    fprintf(fid,'%5.3f',RMS);
    fprintf(fid,'            \n');
    fprintf(fid,'%4.0f%4.0f%4.0f%4.0f\n',[nx;ny;nz;0]);

    % Write out the y/x/z cell spacings
    %=======
    %numbcelly=length(y); %mjc, find number of lines, 7 per line, except last line
    %numbcellx=length(x);
    %numbcellz=length(z);
    %numblinesy=round(numbcelly./7);
    %numblinesx=round(numbcellx./7);
    %numblinesz=round(numbcellz./7);
    %numblinesyend=numbcelly-(numblinesy-1)*7;
    %numblinesxend=numbcellx-(numblinesx-1)*7;
    %numblineszend=numbcellz-(numblinesz-1)*7;
        
    % Write out the x spacings
    for ix=1:length(dx)
        fprintf(fid,'  %9.4E',dx(ix));
        if rem(ix,7)==0 %every 7th model cell write out new line
            fprintf(fid,'\n');
        end
    end
    fprintf(fid,'\n');
    
    % Write out the y spacings
    for iy=1:length(dy) % write in blocks of 7 columns
        fprintf(fid,'  %9.4E',dy(iy)); % 2 spaces only in front
        if rem(iy,7)==0 %every 7th model cell make new line
            fprintf(fid,'\n');
        end
    end
    fprintf(fid,'\n');
    
    % Write out the z spacings
    for iz=1:length(dz)
        fprintf(fid,'  %9.4E',dz(iz));
        if rem(iz,7)==0 %every 7th model cell write out new line
            fprintf(fid,'\n');
        end
    end
    fprintf(fid,'\n');
    
    disp(sprintf('spacings printed')); %#ok<*DSPS>

    for iz=1:nz

        for iy=1:ny

            for ix=nx:-1:1
                fprintf(fid,'  %9.4E',A(ix,iy,iz)); % mjc, add 2 spaces at front of line, also make 4 decimals only
                fprintf(fid,'\n'); % mjc, add this to move to each new line
            end
        end
        disp(sprintf(['rho printed for x,y,z=',num2str(iz)]));
    
        if iz==nz
            disp(sprintf('end line/position reference printed'));
            disp(sprintf('finished printing'));
        end
    end

    fclose(fid);
    close all
end