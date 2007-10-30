function write_adhesion_data(varargin)
% WRITE_STRUCT_DATA     write most the data stored in a given struct to a
%                       set of ascii formated files, using the field
%                       names as the file names
%
%   write_struct_data(S) writes most the field names in struct 's' to ascii
%   formated files, suitable for use in other programs, the fieldnames are
%   used as the file names, all of the files are placed in subfolder
%   'raw_data' of the current working directory
%
%   write_struct_data(S,'dir',d) writes most the field names in struct 's'
%   to ascii formated files, suitable for use in other programs, the
%   fieldnames are used as the file names, all of the files are placed in
%   a subfolder 'raw_data' of the provided directory 'd'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir = 'raw_data';

if (isempty(varargin))
    error('ERROR: write_struct_data - not enough command line options, see ''help write_struct_data''');
elseif (isstruct(varargin{1}))
    data = varargin{1};
else
    error('ERROR: write_struct_data - first parameter was not a struct');
end
    
for i=1:size(varargin,2)
    if (ischar(varargin{i}))
        if (strcmpi(varargin{i},'dir'))
            dir = varargin{i + 1};
        end
    end
end

if (not(exist('data','var')))
    error('ERROR: write_struct_data - can''t seem to find the struct needed');
end

if (not(exist(dir,'dir')))
    mkdir(dir);
end

to_exclude = strvcat('ConvexHull','ConvexImage','Image','FilledImage','PixelIdxList','PixelList','SubarrayIdx');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

field_names = fieldnames(data);

for i = 1:size(field_names,1)
    if(strmatch(field_names(i),to_exclude))
        continue;
    end
    temp = [data.(cell2mat(field_names(i)))];
    save(fullfile(dir,cell2mat(field_names(i))),'temp','-ASCII');
end