function focal_markers = find_focal_adhesion_markers(focal_image,cell_mask);

thresh = adaptive_thresh(focal_image(find(cell_mask)));

focal_markers = im2bw(focal_image,thresh);

focal_markers = imfill(focal_markers,'holes');
focal_markers = bwmorph(focal_markers,'clean');

end