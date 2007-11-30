function edge_binary_image = clean_up_mask_image(input_edge_binary)
%CLEAN_UP_EDGE_IMAGE   isolates the largest enclosed area in a provided
%                       binary image and fills in any holes in the enclosed
%                       area and returns the edge of the filled in area
%
%   EBI = clean_up_mask_image(EB) finds the largest enclosed area in image
%   'EB', isolates that enclosed area, fills in the holes in that area and
%   returns the edge of that filled in area

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;
i_p.FunctionName = 'CLEAN_UP_MASK_IMAGE';

i_p.addRequired('input_edge_binary',@(x) isnumeric(x) || islogical(x));

i_p.parse(input_edge_binary);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
connected_areas = bwlabel(input_edge_binary);
region_sizes = regionprops(connected_areas, 'Area');
max_region_size = max([region_sizes.Area]);
edge_binary_image = ismember(connected_areas, find([region_sizes.Area] == max_region_size));

edge_binary_image = imfill(edge_binary_image,'holes');

edge_binary_image = bwperim(edge_binary_image);

end