function adhesion_props = collect_adhesion_properties(varargin)

%%Error Checking and Setup
if (size(varargin,2) < 1)
    error('ERROR: collect_adhesion_properties - no parameters specified/must specify struct');
elseif (isstruct(varargin{1}))
    image_data = varargin{1};
    if (isfield(image_data,'labeled_adhesions'))
        labeled_adhesions = image_data.labeled_adhesions;
    else
        error('ERROR: collect_adhesion_properties - missing adhesion labels');
    end

    if (isfield(image_data,'original_focal_image'))
        original_focal_image = image_data.original_focal_image;
    else
        error('ERROR: collect_adhesion_properties - missing original focal image');
    end
    
    if (isfield(image_data,'adhesion_props'))
        adhesion_props = image_data.adhesion_props;
    else
        adhesion_props = regionprops(labeled_adhesions,'all');
    end
    
    if (isfield(image_data,'cell_mask'))
        cell_mask = image_data.cell_mask;
    else
        error('ERROR: collect_adhesion_properties - missing original cell mask');
    end
    
else
    error('ERROR: collect_adhesion_properties - first parameter is not a struct');
end

%%Main Program
dists = bwdist(~cell_mask);

for i=1:max(labeled_adhesions(:))
    adhesion_props(i).Average_adhesion_signal = mean(original_focal_image(find(labeled_adhesions == i)));
    adhesion_props(i).Variance_adhesion_signal = var(original_focal_image(find(labeled_adhesions == i)));
    centroid_pos = round(adhesion_props(i).Centroid);
    adhesion_props(i).Centroid_dist_from_edge = dists(centroid_pos(1),centroid_pos(2));
end