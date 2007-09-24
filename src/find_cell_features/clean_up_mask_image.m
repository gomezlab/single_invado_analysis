function edge_binary_image = clean_up_edge_image(input_edge_binary)
%CLEAN_UP_EDGE_IMAGE   isolates the largest enclosed area fills in any
%                      resulting holes 
%
%
connected_areas = bwlabeln(input_edge_binary,8);
region_sizes = regionprops(connected_areas, 'Area');
max_region_size = max([region_sizes.Area]);
edge_binary_image = ismember(connected_areas, find([region_sizes.Area] == max_region_size));

edge_binary_image = imfill(edge_binary_image,'holes');

edge_binary_image = bwperim(edge_binary_image);

end