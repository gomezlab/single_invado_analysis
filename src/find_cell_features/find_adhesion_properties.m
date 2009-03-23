function find_adhesion_properties(focal_file,adhesions_file,varargin)
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
%       -cell_mask: file which contains the cell mask, defaults to not
%        present
%       -debug: set to 1 to output debugging information, defaults to 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'FIND_ADHESION_PROPERTIES';

i_p.addRequired('focal_file',@(x)exist(x,'file') == 2);
i_p.addRequired('adhesions_file',@(x)exist(x,'file') == 2);

i_p.parse(focal_file, adhesions_file);

i_p.addParamValue('output_dir', fileparts(focal_file), @(x)exist(x,'dir')==7);
i_p.addOptional('cell_mask',0,@(x)exist(x,'file') == 2);
i_p.addOptional('protrusion_file',0,@(x)exist(x,'file') == 2);
i_p.addOptional('debug',0,@(x)x == 1 || x == 0);

i_p.parse(focal_file, adhesions_file, varargin{:});

%read in the cell mask image if defined in parameter set
if (isempty(strmatch('cell_mask',i_p.UsingDefaults)))
    cell_mask = imread(i_p.Results.cell_mask);
end

%read in and normalize the input focal adhesion image
focal_image  = imread(focal_file);
scale_factor = double(intmax(class(focal_image)));
focal_image  = double(focal_image)/scale_factor;

%read in the labeled adhesions
adhesions = imread(adhesions_file);

%check if protrusion_file is specified, read it in if specified
if (isempty(strmatch('protrusion_file',i_p.UsingDefaults)))
    load(i_p.Results.protrusion_file);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (exist('cell_mask','var'))
    if (exist('protrusion','var'))
        adhesion_properties = collect_adhesion_properties(focal_image,adhesions,'cell_mask',cell_mask,'protrusion_data',protrusion,'debug',i_p.Results.debug);
    else
        adhesion_properties = collect_adhesion_properties(focal_image,adhesions,'cell_mask',cell_mask,'debug',i_p.Results.debug);
    end
else
    if (exist('protrusion_matrix','var'))
        adhesion_properties = collect_adhesion_properties(focal_image,adhesions,'protrusion_data',protrusion,'debug',i_p.Results.debug);
    else
        adhesion_properties = collect_adhesion_properties(focal_image,adhesions,'debug',i_p.Results.debug);
    end
end
if (i_p.Results.debug), disp('Done with gathering properties'); end

%write the results to files
write_adhesion_data(adhesion_properties,'out_dir',fullfile(i_p.Results.output_dir,'raw_data'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function adhesion_props = collect_adhesion_properties(orig_I,labeled_adhesions,varargin)
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

i_p.addParamValue('cell_mask',0,@(x)isnumeric(x) || islogical(x));
i_p.addParamValue('background_border_size',5,@(x)isnumeric(x));
i_p.addParamValue('protrusion_data',0,@(x)iscell(x));
i_p.addOptional('debug',0,@(x)x == 1 || x == 0);

i_p.parse(labeled_adhesions,orig_I,varargin{:});

%read in the cell mask image if defined in parameter set
if (isempty(strmatch('cell_mask',i_p.UsingDefaults)))
    cell_mask = i_p.Results.cell_mask;
end

if (isempty(strmatch('protrusion_data',i_p.UsingDefaults)))
    protrusion_data = i_p.Results.protrusion_data;
end

adhesion_props = regionprops(labeled_adhesions,'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:max(labeled_adhesions(:))
    adhesion_props(i).Average_adhesion_signal = mean(orig_I(labeled_adhesions == i));
    adhesion_props(i).Variance_adhesion_signal = var(orig_I(labeled_adhesions == i));
    adhesion_props(i).Max_adhesion_signal = max(orig_I(labeled_adhesions == i));
    adhesion_props(i).Min_adhesion_signal = min(orig_I(labeled_adhesions == i));
    
    this_ad = labeled_adhesions;
    this_ad(labeled_adhesions ~= i) = 0;
    this_ad = logical(this_ad);
    background_region = logical(imdilate(this_ad,strel('disk',i_p.Results.background_border_size + 1,0)));
    background_region = and(background_region,not(labeled_adhesions));
    if (exist('cell_mask','var'))
        background_region = and(background_region,cell_mask);
    end
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

if (exist('protrusion_data','var'))
    all_data = [];
    for i=1:size(protrusion_data,2)
        all_data = [all_data, sqrt(protrusion_data{i}(:,3).^2 + protrusion_data{i}(:,4).^2)']; %#ok<AGROW>
    end
    median_velo = median(all_data); %#ok<NASGU>
    
    for i=1:size(protrusion_data,2)
        protrusion_matrix = protrusion_data{i};
        for j=1:max(labeled_adhesions(:))
            dists = sqrt((protrusion_matrix(:,1) - adhesion_props(j).Centroid(1)).^2 + (protrusion_matrix(:,2) - adhesion_props(j).Centroid(2)).^2);
            sorted_dists = sort(dists);
            best_line_nums = find(dists <= sorted_dists(5), 5,'first');
            
            edge_projection = [];
            edge_speed = [];
            for k=1:length(best_line_nums)
                this_line_num = best_line_nums(k);
                
                adhesion_to_edge = [protrusion_matrix(this_line_num,1) - adhesion_props(j).Centroid(1), protrusion_matrix(this_line_num,2) - adhesion_props(j).Centroid(2)];
                adhesion_to_edge = adhesion_to_edge / sqrt(adhesion_to_edge(1)^2 + adhesion_to_edge(2)^2);
                edge_vector = protrusion_matrix(this_line_num,3:4);
                if (sqrt(edge_vector(1)^2 + edge_vector(2)^2) > 10) 
                    edge_vector = (edge_vector / sqrt(edge_vector(1)^2 + edge_vector(2)^2))  * 10;
                end
                
                edge_projection(k) = sqrt(sum(edge_vector.^2))*(dot(edge_vector,adhesion_to_edge)/(sqrt(sum(edge_vector.^2)) * sqrt(sum(adhesion_to_edge.^2)))); %#ok<AGROW>
                edge_speed(k) = sqrt(sum(edge_vector.^2)); %#ok<AGROW>
            end
            adhesion_props(j).Edge_projection(i,1) = mean(edge_projection);
            adhesion_props(j).Edge_speed(i,1) = mean(edge_speed);
        end
    end
end

if (exist('cell_mask','var'))
    [dists, indices] = bwdist(~cell_mask);
    cell_centroid = regionprops(bwlabel(cell_mask),'centroid');
    cell_centroid = cell_centroid.Centroid;
    
    for i=1:max(labeled_adhesions(:))
        centroid_pos = round(adhesion_props(i).Centroid);
        centroid_unrounded = adhesion_props(i).Centroid;
        if(size(centroid_pos,1) == 0)
            warning('MATLAB:noCentroidFound','collect_adhesion_properties - centroid not found');
            adhesion_props(i).Centroid_dist_from_edge = NaN;
        else
            adhesion_props(i).Centroid_dist_from_edge = dists(centroid_pos(2),centroid_pos(1));
            [cep_x,cep_y] = ind2sub(size(cell_mask), indices(centroid_pos(2), centroid_pos(1)));
            adhesion_props(i).Closest_edge_pixel = [cep_x,cep_y];
            
            adhesion_props(i).Centroid_dist_from_center = sqrt((cell_centroid(1) - centroid_unrounded(1))^2 + (cell_centroid(2) - centroid_unrounded(2))^2);
            adhesion_props(i).Angle_to_center = acos((centroid_unrounded(1) - cell_centroid(1))/adhesion_props(i).Centroid_dist_from_center);
            assert(adhesion_props(i).Angle_to_center >= 0 && adhesion_props(i).Angle_to_center <= pi, 'Error: angle to center out of range: %d',adhesion_props(i).Angle_to_center);
            if (centroid_unrounded(2) - cell_centroid(2) < 0)
                if (centroid_unrounded(1) - cell_centroid(1) < 0)
                    assert(adhesion_props(i).Angle_to_center >= pi/2 && adhesion_props(i).Angle_to_center <= pi)
                    adhesion_props(i).Angle_to_center = 2*pi - adhesion_props(i).Angle_to_center;
                elseif (centroid_unrounded(1) - cell_centroid(1) >= 0)
                    assert(adhesion_props(i).Angle_to_center >= 0 && adhesion_props(i).Angle_to_center <= pi/2)
                    adhesion_props(i).Angle_to_center = 2*pi - adhesion_props(i).Angle_to_center;
                end
            end
        end
    end
    
    [border_row,border_col] = ind2sub(size(cell_mask),find(bwperim(cell_mask)));
    adhesion_props(1).Border_pix = [border_col,border_row];
    
    adhesion_props(1).Cell_size = sum(cell_mask(:));
end

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
    'PixelIdxList', 'PixelList', 'SubarrayIdx', 'Border_pix', ...
    'Edge_speed', 'Edge_projection'};

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

print_strings = struct('PixelIdxList','%0.f','Edge_speed','%f','Edge_projection','%f');

for i = 1:size(field_names,1)
    if (not(isempty(strmatch(field_names(i),'PixelIdxList'))) || ...
        not(isempty(strmatch(field_names(i),'Edge_speed')))   || ...
        not(isempty(strmatch(field_names(i),'Edge_projection'))))
        num_ad = size(S,1);
        
        file_out = fullfile(out_dir,[cell2mat(field_names(i)),'.csv']);
        file_handle = fopen(file_out,'wt');
        for j = 1:num_ad
            data = S(j).(field_names{i})';
            for k = 1:size(data,2)
                if (k < size(data,2))
                    fprintf(file_handle,[print_strings.(field_names{i}),','],data(k));
                else
                    fprintf(file_handle,[print_strings.(field_names{i}),'\n'],data(k));
                end
            end
        end
        fclose(file_handle);
    elseif (strmatch(field_names(i),'Border_pix'))
        csvwrite(fullfile(out_dir,[cell2mat(field_names(i)),'.csv']),[S.(cell2mat(field_names(i)))]);
    end
end

% function convert_struct_data_to_mat(S,varargin)

