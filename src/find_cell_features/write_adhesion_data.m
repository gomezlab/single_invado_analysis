function write_adhesion_data(S,varargin)
% WRITE_STRUCT_DATA     write most the data stored in a given struct to a
%                       set of ascii formated files, using the field
%                       names as the file names
%
%   write_struct_data(S) writes most the field names in struct 'S' to ascii
%   formated files, suitable for use in other programs, the fieldnames are
%   used as the file names, all of the files are placed in subfolder
%   'raw_data' of the current working directory
%
%   write_struct_data(S,'out_dir',d) writes most the field names in struct
%   'S' to ascii formated files, suitable for use in other programs, the
%   fieldnames are used as the file names, all of the files are placed in
%   the provided directory 'd'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'WRITE_ADHESION_DATA';

i_p.addRequired('S',@isstruct);
i_p.addParamValue('out_dir','raw_data',@ischar);

i_p.parse(S,varargin{:});
out_dir = i_p.Results.out_dir;

if (not(exist(out_dir,'dir')))
    mkdir(out_dir);
end

to_exclude = strvcat('ConvexHull','ConvexImage','Image','FilledImage','PixelIdxList','PixelList','SubarrayIdx');
special_cases = strvcat('PixelIdxList');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

field_names = fieldnames(S);

for i = 1:size(field_names,1)
    if(strmatch(field_names(i),to_exclude))
        continue;
    end
    temp = [S.(cell2mat(field_names(i)))];
    save(fullfile(out_dir,cell2mat(field_names(i))),'temp','-ASCII');
end

for i = 1:size(field_names,1)
    if (strmatch(field_names(i),special_cases))
        if (strmatch(field_names(i),'PixelIdxList'))
            num_ad = size(S,1);
            if (not(exist(fullfile(out_dir,'PixelIdxList'),'dir')))
                mkdir(fullfile(out_dir,'PixelIdxList'));
            end
            
            for j = 1:num_ad
                temp = [S(j).PixelIdxList];
                save(fullfile(out_dir,'PixelIdxList',num2str(j)),'temp','-ASCII');
            end
        end
    end
end