function edge_binary_image = clean_up_edge_image(input_edge_binary)

connected_areas = bwlabeln(input_edge_binary,8);
region_sizes = regionprops(connected_areas, 'Area');
max_region_size = max([region_sizes.Area]);
edge_binary_image = ismember(connected_areas, find([region_sizes.Area] == max_region_size));

end