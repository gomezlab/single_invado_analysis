function [focal_edge_highlights] = draw_centroid_dots(focal_edge_highlights,cell_mask, watershed_labels)

full_cell_centroid = regionprops(bwlabel(cell_mask),'centroid');
full_cell_centroid = [full_cell_centroid.Centroid];
full_cell_centroid = round(full_cell_centroid);
focal_edge_highlights(full_cell_centroid(2)-5:full_cell_centroid(2)+5,full_cell_centroid(1)-5:full_cell_centroid(1)+5,1) = 1;

centroid_stats = regionprops(watershed_labels,'centroid');
centroid_stats = [centroid_stats.Centroid];

dists = bwdist(~cell_mask);
max_dist = max(dists(:));

for k = 1:2:size(centroid_stats,2)
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