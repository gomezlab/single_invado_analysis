function adhesion_props = collect_adhesion_properties(ad_I,orig_I,varargin)
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

i_p.addRequired('ad_I',@(x)isnumeric(x) || islogical(x));
i_p.addRequired('orig_I',@isnumeric);

i_p.addParamValue('cell_mask',0,@(x)isnumeric(x) || islogical(x));
i_p.addParamValue('background_border_size',5,@(x)isnumeric(x));

i_p.parse(ad_I,orig_I,varargin{:});

%read in the cell mask image if defined in parameter set
if (not(i_p.Results.cell_mask))
else
    cell_mask = i_p.Results.cell_mask;
end

%labeled_adhesions = bwlabel(ad_I,4);
labeled_adhesions = ad_I;
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
    background_region = imdilate(this_ad,strel('square',i_p.Results.background_border_size*2 + 1));
    background_region = and(background_region,not(labeled_adhesions));
    if (exist('cell_mask','var'))
        background_region = and(background_region,cell_mask);
    end
    
    adhesion_props(i).Background_adhesion_signal = mean(orig_I(background_region));
    adhesion_props(i).Background_size = sum(background_region(:));
end

if (exist('cell_mask','var'))
    dists = bwdist(~cell_mask);
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

            adhesion_props(i).Centroid_dist_from_center = sqrt((cell_centroid(1) - centroid_unrounded(1))^2 + (cell_centroid(2) - centroid_unrounded(2))^2);
            adhesion_props(i).Angle_to_center = acos((centroid_unrounded(2) - cell_centroid(2))/adhesion_props(i).Centroid_dist_from_center);
            if (centroid_pos(2) - cell_centroid(2) < 0)
                adhesion_props(i).Angle_to_center = adhesion_props(i).Angle_to_center + pi;
            end
        end

    end
    
    adhesion_props(1).Cell_size = sum(cell_mask(:));
end

end