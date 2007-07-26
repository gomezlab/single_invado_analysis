function focal_edge_highlights = find_each_watershed_adhesion(focal_image, focal_edge_highlights, watershed_labels)

focal_adhesion_sections = zeros(size(focal_image,1),size(focal_image,2));

for i = 1:max(watershed_labels(:))
    thresh = adaptive_thresh(focal_image(find(watershed_labels == i)));

    temp = focal_image;
    temp(find(watershed_labels ~= i)) = 0;

    temp = im2bw(temp,thresh);

    focal_adhesion_sections = temp | focal_adhesion_sections;

end

temp = focal_edge_highlights(:,:,1);

temp(find(focal_adhesion_sections)) = 1;

focal_edge_highlights(:,:,1) = temp;

end