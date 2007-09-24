function adhesion_props = collect_adhesion_properties(varargin)
% COLLECT_ADHESION_PROPERTIES    using the identified adhesions, various
%                                properties are collected concerning the
%                                morphology and physical properties of the
%                                adhesions
%
%   collect_adhesion_properties(I_struct) collects the properties of the
%   adhesions identified in the binary image 'I_struct.adhesions', using
%   the cell mask in 'I_struct.cell_mask' and the original focal image data
%   in 'I_struct.original_image'
%
%   collect_adhesion_properties(ad_I,c_m,orig_I) collects the properties of
%   the adhesions identified in the binary image 'ad_I', using the cell
%   mask in 'c_m' and the original focal image data in 'orig_I'
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
if (isempty(varargin))
    error('ERROR: collect_adhesion_properties - no parameters specified see ''help collect_adhesion_properies''');
elseif (isstruct(varargin{1}))
    image_data = varargin{1};

    if (not(isfield(image_data,'adhesions')))
        error('ERROR: collect_adhesion_properties - missing binary adhesion image, looked in ''adhesions''');
    elseif (not(isfield(image_data,'cell_mask')))
        error('ERROR: collect_adhesion_properties - missing cell mask in struct, looked in ''cell_mask''');
     elseif (not(isfield(image_data,'original_image')))
        error('ERROR: collect_adhesion_properties - missing original focal image, looked in ''original_image''');
    end

    cell_mask = image_data.cell_mask;
    labeled_adhesions = bwlabel(image_data.adhesions);
    original_image = image_data.original_image;

    if (not(isfield(image_data,'adhesion_props')))
        adhesion_props = regionprops(labeled_adhesions,'all');
    else
        adhesion_props = image_data.adhesion_props;
    end

else
    if (lenght(varargin) ~= 3)
        error('ERROR: collect_adhesion_properties - three parameters are required when not using a struct');
    else
        if (not(isnumeric(varargin{1}) || islogical(varargin{1})))
            error('ERROR: collect_adhesion_properties - First parameter should be the binary adhesion image, it appears to not be numeric or logical');
        elseif (not(isnumeric(varargin{2}) || islogical(varargin{2})))
            error('ERROR: collect_adhesion_properties - Second parameter should be the binary cell mask image, it appears to not be numeric or logical');
        elseif (not(isnumeric(varargin{3})))
            error('ERROR: collect_adhesion_properties - Third parameter should be the binary original focal adhesion image, it appears to not be numeric');
        end

        cell_mask = varargin{1};
        labeled_adhesions = bwlabel(varargin{2});
        original_image = varargin{3};
        adhesion_props = regionprops(labeled_adhesions,'all');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dists = bwdist(~cell_mask);
cell_centroid = regionprops(bwlabel(cell_mask),'centroid');
cell_centroid = cell_centroid.Centroid;

for i=1:max(labeled_adhesions(:))
    adhesion_props(i).Average_adhesion_signal = mean(original_image(find(labeled_adhesions == i)));
    adhesion_props(i).Variance_adhesion_signal = var(original_image(find(labeled_adhesions == i)));

    centroid_pos = round(adhesion_props(i).Centroid);
    if(size(centroid_pos,1) == 0)
        warning('collect_adhesion_properties - centroid not found');
        adhesion_props(i).Centroid_dist_from_edge = NaN;
    else
        adhesion_props(i).Centroid_dist_from_edge = dists(centroid_pos(2),centroid_pos(1));
        hypo = sqrt((cell_centroid(1) - centroid_pos(1))^2 + (cell_centroid(2) - centroid_pos(2))^2);
        adhesion_props(i).Angle_to_center = acos((centroid_pos(2) - cell_centroid(2))/hypo);
        if (centroid_pos(2) - cell_centroid(2) < 0)
            adhesion_props(i).Angle_to_center = adhesion_props(i).Angle_to_center + pi;
        end
    end
end