function [focal_edge_highlights] = draw_centroid_dots(image_data)

if (not(isfield(image_data,'focal_edge_highlights')))
    focal_edge_highlights = cat(3, image_data.focal_image,image_data.focal_image,image_data.focal_image);
else
    focal_edge_highlights = image_data.focal_edge_highlights;
end

if (not(isfield(image_data,'cell_mask')))
    error('Missing cell mask data, sorry');
end

if (not(isfield(image_data,'watershed_labels')))
    error('Missing watershed labels data, sorry');
end

full_cell_centroid = regionprops(bwlabel(image_data.cell_mask),'centroid');
full_cell_centroid = [full_cell_centroid.Centroid];
full_cell_centroid = round(full_cell_centroid);
focal_edge_highlights(full_cell_centroid(2)-5:full_cell_centroid(2)+5,full_cell_centroid(1)-5:full_cell_centroid(1)+5,1) = 1;

centroid_stats = regionprops(image_data.watershed_labels,'centroid');
centroid_stats = [centroid_stats.Centroid];

dists = bwdist(~image_data.cell_mask);
max_dist = max(dists(:));

for k = 1:2:size(centroid_stats,2)
    if (isnan(centroid_stats(k)))
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