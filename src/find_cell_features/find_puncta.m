function find_puncta(exp_dir,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.addParamValue('min_puncta_size',1,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('max_puncta_size',Inf,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('max_ratio',Inf,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('filter_size',11,@(x)isnumeric(x) && x > 1);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('..'))
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
filter_thresh = csvread(fullfile(base_dir,image_dirs(1).name,filenames.puncta_threshold));

for i_num = 1:size(image_dirs,1)
    this_image_directory = fullfile(base_dir,image_dirs(i_num).name);
    
    puncta_image = double(imread(fullfile(this_image_directory,filenames.puncta)));
    cell_mask = logical(imread(fullfile(this_image_directory,filenames.cell_mask)));
    
    I_filt = fspecial('disk',i_p.Results.filter_size);
    blurred_image = imfilter(puncta_image,I_filt,'same',mean(puncta_image(:)));
    high_passed_image = puncta_image - blurred_image;
    
    threshed_image = find_threshed_image(high_passed_image,filter_thresh);
    
    %identify and remove objects on the immediate edge of the image
    threshed_image = remove_edge_objects(threshed_image);
    
    puncta = bwlabel(threshed_image,8);
    
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
        filled_ad = imfill(this_puncta,'holes');
        
        puncta(logical(filled_ad)) = this_num;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Size and Eccentricity Filtering
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
    
    %     if (not(any(strcmp('max_solidity',i_p.UsingDefaults))))
    %         props = regionprops(puncta,'Solidity');
    %         puncta = puncta .* ismember(puncta, find([props.Solidity] <= i_p.Results.max_eccentricity));
    %     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Renumber objects to be sequential
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ad_nums = unique(puncta);
    assert(ad_nums(1) == 0, 'Background pixels not found after building puncta label matrix')
    for i = 2:length(ad_nums)
        puncta(puncta == ad_nums(i)) = i - 1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Build puncta perimeters image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    puncta_perim = zeros(size(puncta));
    for i = 1:max(puncta(:))
        assert(any(any(puncta == i)), 'Error: can''t find ad number %d', i);
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
    
    scaled_image = (puncta_image - min_max(1))/(range(min_max));
    scaled_image(scaled_image < 0) = 0;
    scaled_image(scaled_image > 1) = 1;
    highlighted_image = create_highlighted_image(scaled_image, im2bw(puncta_perim), ...
        'color_map',[1 0 0]);
    if (exist('cell_mask','var'))
        highlighted_image = create_highlighted_image(highlighted_image, bwperim(cell_mask), ...
            'color_map',[0 1 0],'mix_percent',0.5);
    end
    if (size(highlighted_image,1) > 800)
        highlighted_image = imresize(highlighted_image,[800,NaN]);
    end
    imwrite(highlighted_image, fullfile(this_image_directory, filenames.objects_highlight));
    
    if (mod(i_num,10)==0)
        disp(['Done with ',num2str(i_num),'/',num2str(size(image_dirs,1))])
    end
end

toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function threshed_image = find_threshed_image(high_passed_image, filter_thresh)

if (length(filter_thresh) == 1)
    threshed_image = high_passed_image >= filter_thresh;
else
    high_threshed_image = high_passed_image >= filter_thresh(2);
    high_threshed_image = remove_edge_objects(high_threshed_image);
    
    low_threshed_image = high_passed_image >= filter_thresh(1);
    low_thresh_bwlabel = bwlabel(low_threshed_image,4);
    
    overlap_labels = unique(low_thresh_bwlabel.*high_threshed_image);
    if (overlap_labels(1) == 0)
        overlap_labels = overlap_labels(2:end);
    end
    
    threshed_image = ismember(low_thresh_bwlabel,overlap_labels);
end