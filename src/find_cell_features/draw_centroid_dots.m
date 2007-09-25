function [focal_edge_highlights] = draw_centroid_dots(image_data)
% DRAW_CENTROID_DOTS    using a given labeled watershed image, dots are
%                       overlayed on the original focal image, the color of
%                       the dots depends on the distance of the centroid of
%                       the labeled adhesion from the edge of the image
%
%   draw_centroid_dots(I_struct) using the image in
%   'I_struct.focal_edge_highlights' (if not present, the image in
%   'I_struct.original_focal_image' is used as a base) the centroids of
%   each region in 'I_struct.watershed labels' (if not present, they are
%   calculated from 'I_struct.watershed_edges'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (not(isfield(image_data,'focal_edge_highlights')))
    focal_edge_highlights = cat(3, image_data.original_focal_image,image_data.original_focal_image,image_data.original_focal_image);
else
    focal_edge_highlights = image_data.focal_edge_highlights;
end

if (not(isfield(image_data,'cell_mask')))
    error('ERROR: draw_centroid_dots - Missing cell mask data, sorry');
else
    cell_mask = image_data.cell_mask;
end

if (not(isfield(image_data,'watershed_labels')))
    if (not(isfield(image_data,'watershed_edges')))
        error('ERROR: draw_centroid_dots - Missing watershed labels and watershed edges data, sorry');
    else
        watershed_labels = bwlabel(~watershed_edges.*cell_mask);
    end
else
    watershed_labels = image_data.watereshed_labels;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
full_cell_centroid = regionprops(bwlabel(image_data.cell_mask),'Centroid');
full_cell_centroid = [full_cell_centroid.Centroid];
full_cell_centroid = round(full_cell_centroid);
focal_edge_highlights(full_cell_centroid(2)-5:full_cell_centroid(2)+5,full_cell_centroid(1)-5:full_cell_centroid(1)+5,1) = 1;

centroid_stats = regionprops(watershed_labels,'Centroid');
centroid_stats = [centroid_stats.Centroid];

dists = bwdist(~image_data.cell_mask);
max_dist = max(dists(:));

for k = 1:2:size(centroid_stats,2)
    if (isnan(centroid_stats(k)))
        warning('draw_centroid_dots - NaN encountered in centroid value');
        continue;
    end

    middle_point_col = round(centroid_stats(k));
    middle_point_row = round(centroid_stats(k+1));
    col_range = middle_point_col-1:middle_point_col+1;
    row_range = middle_point_row-1:middle_point_row+1;

    this_dist = dists(middle_point_row,middle_point_col)/max_dist;

    color = ones(3,3,3);
    color(:,:,1) = 1-this_dist;
    color(:,:,2) = zeros(3,3);
    color(:,:,3) = this_dist;

    focal_edge_highlights(row_range,col_range,1:3) = color;
end

end