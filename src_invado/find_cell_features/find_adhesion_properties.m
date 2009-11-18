function find_adhesion_properties(focal_file,adhesions_file,gel_binary_file,varargin)
% FIND_ADHESION_PROPERTIES    deteremines and outputs the quantitative
%                             properties associated with the adhesions
%                             located in prior steps
%
%   find_adhesion_properties(ff,af,OPTIONS) determines the quantitative
%   properites of the adhesions identified in the file 'af', using
%   information in the focal image file 'ff', the properties are written
%   to a set of csv files located a folder named 'raw_data' in the same
%   directory as the focal image file 'ff'
%
%   Options:
%
%       -debug: set to 1 to output debugging information, defaults to 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'FIND_ADHESION_PROPERTIES';

i_p.addRequired('focal_file',@(x)exist(x,'file') == 2);
i_p.addRequired('adhesions_file',@(x)exist(x,'file') == 2);
i_p.addRequired('gel_binary_file',@(x)exist(x,'file') == 2);

i_p.parse(focal_file, adhesions_file,gel_binary_file);

i_p.addParamValue('output_dir', fileparts(focal_file), @(x)exist(x,'dir')==7);
i_p.addOptional('debug',0,@(x)x == 1 || x == 0);

i_p.parse(focal_file, adhesions_file, gel_binary_file, varargin{:});

%read in and normalize the input focal adhesion image
focal_image  = imread(focal_file);
scale_factor = double(intmax(class(focal_image)));
focal_image  = double(focal_image)/scale_factor;

%read in the labeled adhesions
adhesions = imread(adhesions_file);

%read in the labeled adhesions
gel_binary = logical(imread(gel_binary_file));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
adhesion_properties = collect_adhesion_properties(focal_image,adhesions,gel_binary,'debug',i_p.Results.debug);
if (i_p.Results.debug), disp('Done with gathering properties'); end

%write the results to files
write_adhesion_data(adhesion_properties,'out_dir',fullfile(i_p.Results.output_dir,'raw_data'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function adhesion_props = collect_adhesion_properties(orig_I,labeled_adhesions,gel_binary,varargin)
% COLLECT_ADHESION_PROPERTIES    using the identified adhesions, various
%                                properties are collected concerning the
%                                morphology and physical properties of the
%                                adhesions
%
%   ad_p = collect_adhesion_properties(ad_I,c_m,orig_I) collects the
%   properties of the adhesions identified in the binary image 'ad_I',
%   using the cell mask in 'c_m' and the original focal image data in
%   'orig_I', returning a structure 'ad_p' containing properties
%
%   Properties Collected:
%       -all of the properties collected by regioprops(...,'all')
%       -the distance of each adhesion's centroid from the nearest cell
%        edge
%       -the average and variance of the normalized fluorescence signal
%        within each adhesion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;
i_p.FunctionName = 'COLLECT_ADHESION_PROPERTIES';

i_p.addRequired('orig_I',@isnumeric);
i_p.addRequired('labeled_adhesions',@(x)isnumeric(x));
i_p.addRequired('gel_binary',@(x)islogical(x));

i_p.addParamValue('background_border_size',5,@(x)isnumeric(x));
i_p.addOptional('debug',0,@(x)x == 1 || x == 0);

i_p.parse(labeled_adhesions,orig_I,gel_binary,varargin{:});

adhesion_props = regionprops(labeled_adhesions,'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Properites Always Extracted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:max(labeled_adhesions(:))
    adhesion_props(i).Average_adhesion_signal = mean(orig_I(labeled_adhesions == i));
    adhesion_props(i).Variance_adhesion_signal = var(orig_I(labeled_adhesions == i));
    adhesion_props(i).Max_adhesion_signal = max(orig_I(labeled_adhesions == i));
    adhesion_props(i).Min_adhesion_signal = min(orig_I(labeled_adhesions == i));
    
    this_ad = labeled_adhesions;
    this_ad(labeled_adhesions ~= i) = 0;
    this_ad = logical(this_ad);
    background_region = logical(imdilate(this_ad,strel('disk',i_p.Results.background_border_size,0)));
    background_region = and(background_region,not(labeled_adhesions));
    
    adhesion_props(i).Degrade_overlap = any(any(this_ad .* gel_binary));
    adhesion_props(i).Degrade_overlap_percent = sum(sum(this_ad .* gel_binary))/adhesion_props(i).Area;
    assert(adhesion_props(i).Degrade_overlap_percent >= 0 && adhesion_props(i).Degrade_overlap_percent <= 1);
    
    adhesion_props(i).Background_adhesion_signal = mean(orig_I(background_region));
    adhesion_props(i).Background_area = sum(background_region(:));
    adhesion_props(i).Background_corrected_signal = adhesion_props(i).Average_adhesion_signal - adhesion_props(i).Background_adhesion_signal;
    
    shrunk_region = logical(imerode(this_ad,strel('disk',1,0)));
    if (sum(shrunk_region(:)) == 0), shrunk_region = this_ad; end
    adhesion_props(i).Shrunk_area = sum(shrunk_region(:));
    adhesion_props(i).Shrunk_adhesion_signal = mean(orig_I(shrunk_region));
    adhesion_props(i).Shrunk_corrected_signal = adhesion_props(i).Shrunk_adhesion_signal - adhesion_props(i).Background_adhesion_signal;
    
    if (mod(i,10) == 0 && i_p.Results.debug), disp(['Finished Ad: ',num2str(i), '/', num2str(max(labeled_adhesions(:)))]); end
end

adhesion_mask = im2bw(labeled_adhesions,0);
adhesion_props(1).Adhesion_mean_intensity = sum(sum(orig_I(adhesion_mask)))/sum(sum(adhesion_mask));

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

to_exclude = {'ConvexHull','ConvexImage','Image','FilledImage', ...
    'PixelList', 'SubarrayIdx', 'Border_pix', 'Extrema'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

field_names = fieldnames(S);

ad_props_cell = struct2cell(S);

print_strings = struct('PixelIdxList','%0.f');

for i = 1:size(field_names,1)
    
    if(strmatch(field_names(i),to_exclude))
        continue;
    end
    
    format_string = '%f';
    if(strmatch(field_names(i),fieldnames(print_strings)))
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
