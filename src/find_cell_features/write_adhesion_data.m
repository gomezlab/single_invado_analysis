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

to_exclude = {'ConvexHull','ConvexImage','Image','FilledImage','PixelIdxList','PixelList','SubarrayIdx', 'Border_pix'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

field_names = fieldnames(S);

for i = 1:size(field_names,1)

    if(strmatch(field_names(i),to_exclude))
        continue;
    end

    data = [S.(cell2mat(field_names(i)))];

    file_out = fullfile(out_dir,[cell2mat(field_names(i)),'.csv']);
    file_handle = fopen(file_out,'wt');
    for j = 1:size(data,2)
        if (j < size(data,2))
            fprintf(file_handle,'%f,',data(j));
        else
            fprintf(file_handle,'%f',data(j));
        end
    end
    fclose(file_handle);
end

for i = 1:size(field_names,1)
    if (strmatch(field_names(i),'PixelIdxList'))
        num_ad = size(S,1);

        file_out = fullfile(out_dir,[cell2mat(field_names(i)),'.csv']);
        file_handle = fopen(file_out,'wt');
        for j = 1:num_ad
            data = (S(j).PixelIdxList');
            for k = 1:size(data,2)
                if (k < size(data,2))
                    fprintf(file_handle,'%0.f,',data(k));
                else
                    fprintf(file_handle,'%0.f\n',data(k));
                end
            end
        end
        fclose(file_handle);
    elseif (strmatch(field_names(i),'Border_pix'))
        csvwrite(fullfile(out_dir,[cell2mat(field_names(i)),'.csv']),[S.(cell2mat(field_names(i)))]);
    end
end
