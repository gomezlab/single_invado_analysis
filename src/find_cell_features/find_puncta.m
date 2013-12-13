function find_puncta(exp_dir,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.addParamValue('stdev_thresh',[1.75,3],@(x)isnumeric(x));
i_p.addParamValue('min_puncta_seed_size',3,@(x)isnumeric(x));
i_p.addParamValue('min_puncta_size',1,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('max_puncta_size',Inf,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('max_ratio',Inf,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('filter_size',11,@(x)isnumeric(x) && x > 1);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'))
addpath(genpath('../visualize_cell_features'))
filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_dir = fullfile(exp_dir,'individual_pictures');

image_dirs = dir(base_dir);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);

min_max = csvread(fullfile(base_dir,image_dirs(1).name,filenames.puncta_range));
thresh_settings = csvread(fullfile(base_dir,image_dirs(1).name,filenames.puncta_threshold));
filter_thresh = thresh_settings(1) + i_p.Results.stdev_thresh*thresh_settings(2);

for i_num = 1:size(image_dirs,1)
    this_image_directory = fullfile(base_dir,image_dirs(i_num).name);
    
    puncta_image = double(imread(fullfile(this_image_directory,filenames.puncta)));
    gel_image = double(imread(fullfile(this_image_directory,filenames.gel)));
    cell_mask = logical(imread(fullfile(this_image_directory,filenames.cell_mask)));
    
    I_filt = fspecial('disk',i_p.Results.filter_size);
    blurred_image = imfilter(puncta_image,I_filt,'replicate');
    high_passed_image = puncta_image - blurred_image;
    
    threshed_image = find_threshed_image(high_passed_image,filter_thresh,...
        i_p.Results.min_puncta_seed_size);
    
    %identify and remove objects on the immediate edge of the image
    threshed_image = remove_edge_objects(threshed_image);
    
    puncta = bwlabel(threshed_image,8);
    imwrite(threshed_image,fullfile(this_image_directory, filenames.objects_above_thresh));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Remove objects outside mask
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (exist('cell_mask','var'))
        inside_cell_mask = unique(puncta(cell_mask));
        inside_cell_mask = inside_cell_mask(inside_cell_mask > 0);
        
        puncta = puncta .* ismember(puncta,inside_cell_mask);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find and fill holes in single puncta
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    props = regionprops(puncta,'Area');
    large_punc_nums = find([props.Area] >= 4);
    for this_num = large_punc_nums
        %first make a binary image of the current object and then run imfill
        %to fill any holes present
        this_puncta = puncta == this_num;
        filled_obj = imfill(this_puncta,'holes');
        
        %dealing with an edge case where one puncta completely surrounds
        %other objects, in that case, make sure those other objects don't
        %get reassigned based on being enclosed
        other_puncta = puncta ~= this_num & puncta ~= 0;
        filled_obj(other_puncta) = 0;
        
        puncta(logical(filled_obj)) = this_num;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Size and Ratio Filtering
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (i_p.Results.min_puncta_size > 1)
        props = regionprops(puncta,'Area');
        puncta = puncta .* ismember(puncta, find([props.Area] >= i_p.Results.min_puncta_size));
    end
    if (not(isinf(i_p.Results.max_puncta_size)))
        props = regionprops(puncta,'Area');
        puncta = puncta .* ismember(puncta, find([props.Area] <= i_p.Results.max_puncta_size));
    end
    
    if (not(isinf(i_p.Results.max_ratio)))
        props = regionprops(puncta,'MajorAxisLength','MinorAxisLength');
        ratio = [props.MajorAxisLength]./[props.MinorAxisLength];
        puncta = puncta .* ismember(puncta, find(ratio <= i_p.Results.max_ratio));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Renumber objects to be sequential
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj_nums = unique(puncta);
    assert(obj_nums(1) == 0, 'Background pixels not found after building puncta label matrix')
    for i = 2:length(obj_nums)
        puncta(puncta == obj_nums(i)) = i - 1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Build puncta perimeters image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    puncta_perim = zeros(size(puncta));
    for i = 1:max(puncta(:))
        assert(any(any(puncta == i)), 'Error: can''t find obj number %d', i);
        this_puncta = zeros(size(puncta));
        this_puncta(puncta == i) = 1;
        puncta_perim(bwperim(this_puncta)) = i;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Write the output files
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    imwrite(double(puncta)/2^16,fullfile(this_image_directory, filenames.objects),'bitdepth',16);
    imwrite(double(puncta_perim)/2^16,fullfile(this_image_directory, filenames.objects_perim),'bitdepth',16);
    imwrite(im2bw(puncta),fullfile(this_image_directory, filenames.objects_binary));
    
    gel_norm = (gel_image - min(gel_image(:)))/range(gel_image(:));
    
    scaled_image = (puncta_image - min_max(1))/(range(min_max));
    scaled_image(scaled_image < 0) = 0;
    scaled_image(scaled_image > 1) = 1;
    highlighted_image = create_highlighted_image(scaled_image, im2bw(puncta_perim), ...
        'color_map',[1 0 0]);
    highlighted_gel = create_highlighted_image(gel_norm, im2bw(puncta_perim), ...
        'color_map',[1 0 0]);
    
    if (exist('cell_mask','var'))
        highlighted_image = create_highlighted_image(highlighted_image, bwperim(cell_mask), ...
            'color_map',[0 1 0],'mix_percent',0.5);
        highlighted_gel = create_highlighted_image(highlighted_gel, bwperim(cell_mask), ...
            'color_map',[0 1 0],'mix_percent',0.5);
    end
    imwrite(highlighted_image, fullfile(this_image_directory, filenames.objects_highlight));
    imwrite(highlighted_gel, fullfile(this_image_directory, 'gel_highlight.png'));
    
    if (mod(i_num,10)==0)
        disp(['Done with ',num2str(i_num),'/',num2str(size(image_dirs,1))])
    end
end

toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function threshed_image = find_threshed_image(high_passed_image, filter_thresh, min_seed_size)

if (length(filter_thresh) == 1)
    threshed_image = high_passed_image >= filter_thresh;
    threshed_image = remove_edge_objects(threshed_image);
else
    high_threshed_image = high_passed_image >= filter_thresh(2);
    high_threshed_image = remove_edge_objects(high_threshed_image);
    high_threshed_image_label = bwlabel(high_threshed_image,4);
    high_threshed_image_label = filter_label_mat_by_size(high_threshed_image_label,min_seed_size);
    high_threshed_image = high_threshed_image_label > 0;
    
    low_threshed_image = high_passed_image >= filter_thresh(1);
    low_thresh_bwlabel = bwlabel(low_threshed_image,4);
    
    overlap_labels = unique(low_thresh_bwlabel.*high_threshed_image);
    overlap_labels = overlap_labels(overlap_labels > 0);
    
    threshed_image = ismember(low_thresh_bwlabel,overlap_labels);
    threshed_image = remove_edge_objects(threshed_image);
end

function filt_image = filter_label_mat_by_size(label_mat,size_mat)

props = regionprops(label_mat,'Area');

if (length(size_mat) == 1)
    filt_image = ismember(label_mat,find([props.Area] >= size_mat(1))).*label_mat;
else
    filt_image = ismember(label_mat,find([props.Area] >= size_mat(1) & [props.Area] <= size_mat(2))).*label_mat;
end