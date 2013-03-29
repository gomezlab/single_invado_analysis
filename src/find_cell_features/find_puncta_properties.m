function find_puncta_properties(exp_dir,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.addOptional('write_to_flat_files',0,@(x)x == 1 | x == 0);
i_p.addOptional('debug',0,@(x)x == 1 | x == 0);

i_p.parse(exp_dir,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'));

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base_dir = fullfile(exp_dir,'individual_pictures');

image_dirs = dir(base_dir);
image_dirs = image_dirs(3:end);

all_cell_props = cell(0);

for i_num = 1:size(image_dirs,1)
    current_data = read_in_file_set(fullfile(base_dir,image_dirs(i_num).name),filenames);
    if (unique(current_data.objects(:)) == 0)
        continue;
    end
    object_properties = collect_object_properties(current_data,'debug',i_p.Results.debug);
    all_cell_props{i_num} = object_properties;
    
    %write the results to files
    if (i_p.Results.write_to_flat_files)
        write_object_data(object_properties,...
            'out_dir',fullfile(base_dir,image_dirs(i_num).name,'raw_data'));
    end
    
    if (mod(i_num,10)==0)
        disp(['Done with ',num2str(i_num),'/',num2str(size(image_dirs,1))])
    end
end

save(fullfile(base_dir, image_dirs(1).name,filenames.cell_props),'all_cell_props');

toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function object_props = collect_object_properties(c_d,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('c_d',@isstruct);

i_p.addOptional('debug',0,@(x)x == 1 || x == 0);

i_p.parse(c_d,varargin{:});

object_props = regionprops(c_d.objects,'Area','Centroid','Eccentricity');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Properites Always Extracted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:max(c_d.objects(:))
    object_props(i).Average_puncta_signal = mean(c_d.puncta_image(c_d.objects == i));
    
    %this bit of code isolates a single object as a logical image and
    %then builds another logical image of the region around the object,
    %excluding certain areas
    this_ad = c_d.objects;
    this_ad(c_d.objects ~= i) = 0;
    this_ad = logical(this_ad);
    
    current_diffs = collect_local_diff_properties(c_d,this_ad);
    
    diff_names = fieldnames(current_diffs);
    for j=1:length(diff_names)
        object_props(i).(diff_names{j}) = current_diffs.(diff_names{j});
    end
        
    if (mod(i,10) == 0 && i_p.Results.debug), disp(['Finished Ad: ',num2str(i), '/', num2str(max(c_d.objects(:)))]); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Properites Extracted If Cell Mask Available
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (isfield(c_d, 'cell_mask'))
    [dists, ~] = bwdist(~c_d.cell_mask);
    
    %Now we search for the pixels which are closest to an edge of the cell
    %mask that is also touching the edge of image. We want to find these
    %pixels because the true closest cell edge may be off the microscope
    %field of view. To be safe, we will set those
    %distance-to-nearest-cell-edge values to NaN.
    black_border_mask = ones(size(c_d.cell_mask,1),size(c_d.cell_mask,2));
    black_border_mask(1,:) = 0; black_border_mask(end,:) = 0;
    black_border_mask(:,1) = 0; black_border_mask(:,end) = 0;
    
    bb_dists_temp = bwdist(~black_border_mask);
        
    dists(bb_dists_temp < dists) = NaN;
    for i=1:max(c_d.objects(:))
        centroid_pos = round(object_props(i).Centroid);
        centroid_unrounded = object_props(i).Centroid; %#ok<NASGU>
        if(size(centroid_pos,1) == 0)
            warning('MATLAB:noCentroidFound','collect_object_properties - centroid not found');
            object_props(i).Centroid_dist_from_edge = NaN;
        else
            object_props(i).Centroid_dist_from_edge = dists(centroid_pos(2),centroid_pos(1));
        end
    end
end

function write_object_data(S,varargin)
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

i_p.addRequired('S',@isstruct);
i_p.addParamValue('out_dir','raw_data',@ischar);

i_p.parse(S,varargin{:});
out_dir = i_p.Results.out_dir;

if (not(exist(out_dir,'dir')))
    mkdir(out_dir);
end

to_exclude = {'ConvexHull','ConvexImage','Image','FilledImage', ...
    'PixelList', 'SubarrayIdx', 'Border_pix', 'Extrema'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

field_names = fieldnames(S);

ad_props_cell = struct2cell(S);

print_strings = struct('PixelIdxList','%0.f',...
    'Orientation','%0.2f','MajorAxisLength','%0.2f',...
    'MinorAxisLength','%0.2f','Area','%d','Centroid_dist_from_edge','%0.2f',...
    'Local_gel_diff','%0.2f');

for i = 1:size(field_names,1)
    if(strcmp(field_names(i),to_exclude))
        continue;
    end
    
    format_string = '%f';
    if(strcmp(field_names(i),fieldnames(print_strings)))
        format_string = print_strings.(field_names{i});
    end
    
    file_out = fullfile(out_dir,[cell2mat(field_names(i)),'.csv']);
    
    data = ad_props_cell(i,:);
    output_CSV_from_cell(data, file_out, 'format', format_string);
end

function output_CSV_from_cell(data, out_file, varargin)
% output_CSV_from_cell    writes a provided cell data structure to a CSV
%                         file
%
%   write_struct_data(D,OF) writes the data in cell 'D' to file 'DF', using
%   fprintf format '%f'
%
%   Options:
%       -format: specify the string passed to fprintf for number conversion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'output_CSV_from_cell';

i_p.addRequired('data',@iscell);
i_p.addRequired('out_file', @(x) exist(fileparts(x),'dir') == 7 );

i_p.addParamValue('format','%f',@ischar);

i_p.parse(data,out_file,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_handle = fopen(i_p.Results.out_file,'wt');
for i = 1:max(size(data))
    for j = 1:max(size(data{i}))
        assert(any(size(data{i}) == 1))
        if (j < max(size(data{i})))
            fprintf(file_handle,[i_p.Results.format,','],data{i}(j));
        else
            fprintf(file_handle,[i_p.Results.format,'\n'],data{i}(j));
        end
    end
end
fclose(file_handle);