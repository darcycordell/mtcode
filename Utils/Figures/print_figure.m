function print_figure(folder, filename, override)
%
% Function which saves a figure. The format of the figure is specified in
% user_defaults.m as png, eps, etc.
%
% Usage: print_figure(folder,filename,override)
%
% Inputs: 
%       folder: specfies the folder to save the file in. If the folder does
%       not exist, the folder is created. You can specify a full path as a
%       string or just a folder name (e.g. 'C:/User/Plots/' or 'Plots'). If
%       you do not specify a full path then the folder is created in the
%       current directory
%
%       filename: A string filename for the figure file
%
%       override: An OPTIONAL input (1 or 0) to override user_defaults and
%       automatically output a PNG file. This is used in NIMS processing to
%       avoid accidentally outputting large numbers of EPS files by
%       accident.
%

if ~exist('override','var')
    override = 0;
end

u = user_defaults;

if strcmp(u.output_figure,'none')
    return
end

if override
    u.output_figure = 'png';
end
    
% output = string(u.output_figure); % put all input file types in a string array but only works in 2016b and later...
output = cellstr(u.output_figure); % this is a cell array 

for ifig = 1:numel(output)

    if strcmp(char(output(ifig)),'jpeg') || strcmp(char(output(ifig)),'JPEG') || strcmp(char(output(ifig)),'jpg') || strcmp(char(output(ifig)),'JPG')
        
        if exist(['./',folder],'dir')~=7
            mkdir(folder);
        end

        printfile=[folder,'/',filename];
        print('-djpeg',[printfile,'.jpeg'])

    elseif strcmp(char(output(ifig)),'eps') || strcmp(char(output(ifig)),'EPS')
        
        if u.plot_topo
            check_menu = menu('WARNING: Saving an EPS file with SRTM topography plotted may be very slow and may crash MATLAB. Are you sure you want to proceed?','Yes','No');
            if check_menu == 2
                disp('WARNING: Saving an EPS file with SRTM topography plotted may be very slow and may crash MATLAB');
                disp('Please check user_defaults and either change output_figure type or plot_topo = false')
                return
            end
        end
        
        if exist(['./',folder],'dir')~=7
            mkdir(folder);
        end

        printfile=[folder,'/',filename];
        print('-depsc','-painters',[printfile,'.eps'])

    elseif strcmp(char(output(ifig)),'png') || strcmp(char(output(ifig)),'PNG')
        
        if exist(['./',folder],'dir')~=7
            mkdir(folder);
        end

        printfile=[folder,'/',filename];
        print('-dpng',[printfile,'.png'])
        
    else
        
        disp('Unrecognized file type in u.output_figure. Please enter .png, .jpg, or .eps')

    end

end


end %END print_figure