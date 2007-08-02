function focal_markers = find_focal_adhesion_markers(image_data);

if (not(isfield(image_data,'focal_image')))
    error('Missing focal image data, sorry');
end
if (not(isfield(image_data,'cell_mask')))
    error('Missing cell mask data, sorry');
end

thresh = adaptive_thresh(image_data.focal_image(find(image_data.cell_mask)));

focal_markers = im2bw(image_data.focal_image,thresh);

focal_markers = imfill(focal_markers,'holes');
focal_markers = bwmorph(focal_markers,'clean');

end