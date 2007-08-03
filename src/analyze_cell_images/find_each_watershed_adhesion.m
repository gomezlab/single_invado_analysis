function [focal_edge_highlights,varargout] = find_watershed_adhesions(focal_image, focal_edge_highlights, watershed_labels, cell_mask)

if (not(isfield(image_data,'focal_edge_highlights')))
    focal_edge_highlights = cat(3, image_data.focal_image,image_data.focal_image,image_data.focal_image);
else
    focal_edge_highlights = image_data.focal_edge_highlights;
end

if (not(isfield(image_data,'focal_image')))
    error('Missing cell mask data, sorry');
end
if (not(isfield(image_data,'cell_mask')))
    error('Missing cell mask data, sorry');
end
if (not(isfield(image_data,'watershed_labels')))
    error('Missing watershed labels data, sorry');
end

focal_adhesion_sections = zeros(size(focal_image,1),size(focal_image,2));
min_thresh = adaptive_thresh(focal_image(find(cell_mask)));

for i = 1:max(watershed_labels(:))
    thresh = adaptive_thresh(focal_image(find(watershed_labels == i)));
    
    if (thresh < min_thresh)
        thresh = min_thresh;
    end
    temp = focal_image;
    temp(find(watershed_labels ~= i)) = 0;

    temp = im2bw(temp,thresh);

    focal_adhesion_sections = temp | focal_adhesion_sections;

end

if (nargout > 1)
    varargout(1) = {focal_adhesion_sections};
end

temp = focal_edge_highlights(:,:,1);

temp(find(focal_adhesion_sections)) = 1;

focal_edge_highlights(:,:,1) = temp;

end