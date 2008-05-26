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

i_p.parse(ad_I,orig_I,varargin{:});

%read in the cell mask image if defined in parameter set
if (not(i_p.Results.cell_mask))
else
    cell_mask = i_p.Results.cell_mask;
end

labeled_adhesions = bwlabel(ad_I,4);
adhesion_props = regionprops(labeled_adhesions,'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Average_adhesion_signal = zeros(1,max(labeled_adhesions(:)));
Variance_adhesion_signal = zeros(1,max(labeled_adhesions(:)));
Max_adhesion_signal = zeros(1,max(labeled_adhesions(:)));
Min_adhesion_signal = zeros(1,max(labeled_adhesions(:)));

for i=1:max(labeled_adhesions(:))
    Average_adhesion_signal(i) = mean(orig_I(labeled_adhesions == i));
    Variance_adhesion_signal(i) = var(orig_I(labeled_adhesions == i));
    Max_adhesion_signal(i) = max(orig_I(labeled_adhesions == i));
    Min_adhesion_signal(i) = min(orig_I(labeled_adhesions == i));
end

adhesion_props(1).Average_adhesion_signal = Average_adhesion_signal;
adhesion_props(1).Average_pixel_adhesion_signal = Average_adhesion_signal*(2^12 - 1);
adhesion_props(1).Variance_adhesion_signal = Variance_adhesion_signal;
adhesion_props(1).Max_adhesion_signal = Max_adhesion_signal;
adhesion_props(1).Min_adhesion_signal = Min_adhesion_signal;

if (exist('cell_mask','var'))
    dists = bwdist(~cell_mask);
    cell_centroid = regionprops(bwlabel(cell_mask),'centroid');
    cell_centroid = cell_centroid.Centroid;

    Centroid_dist_from_edge = zeros(1,max(labeled_adhesions(:)));
    Centroid_dist_from_center = zeros(1,max(labeled_adhesions(:)));
    Angle_to_center = zeros(1,max(labeled_adhesions(:)));

    for i=1:max(labeled_adhesions(:))
        Average_adhesion_signal(i) = mean(orig_I(labeled_adhesions == i));
        Variance_adhesion_signal(i) = var(orig_I(labeled_adhesions == i));

        centroid_pos = round(adhesion_props(i).Centroid);
        centroid_unrounded = adhesion_props(i).Centroid;
        if(size(centroid_pos,1) == 0)
            warning('MATLAB:noCentroidFound','collect_adhesion_properties - centroid not found');
            adhesion_props(i).Centroid_dist_from_edge = NaN;
        else
            Centroid_dist_from_edge(i) = dists(centroid_pos(2),centroid_pos(1));

            Centroid_dist_from_center(i) = sqrt((cell_centroid(1) - centroid_unrounded(1))^2 + (cell_centroid(2) - centroid_unrounded(2))^2);
            Angle_to_center(i) = acos((centroid_unrounded(2) - cell_centroid(2))/Centroid_dist_from_center(i));
            if (centroid_pos(2) - cell_centroid(2) < 0)
                Angle_to_center(i) = Angle_to_center(i) + pi;
            end
        end

    end

    adhesion_props(1).Centroid_dist_from_edge = Centroid_dist_from_edge;
    adhesion_props(1).Centroid_dist_from_center = Centroid_dist_from_center;
    adhesion_props(1).Angle_to_center = Angle_to_center;
    adhesion_props(1).Cell_size = sum(cell_mask(:));
end

end