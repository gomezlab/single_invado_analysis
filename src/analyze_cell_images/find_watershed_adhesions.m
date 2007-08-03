function focal_adhesion_sections = find_watershed_adhesions(image_data)

if (not(isfield(image_data,'focal_image')))
    error('Missing cell mask data, sorry');
end
if (not(isfield(image_data,'cell_mask')))
    error('Missing cell mask data, sorry');
end
if (not(isfield(image_data,'watershed_labels')))
    error('Missing watershed labels data, sorry');
end

focal_adhesion_sections = zeros(size(image_data.focal_image,1),size(image_data.focal_image,2));
min_thresh = adaptive_thresh(image_data.focal_image(find(image_data.cell_mask)));

for i = 1:max(image_data.watershed_labels(:))
    thresh = adaptive_thresh(image_data.focal_image(find(image_data.watershed_labels == i)));
    
    if (thresh < min_thresh)
        thresh = min_thresh;
    end
    temp = image_data.focal_image;
    temp(find(image_data.watershed_labels ~= i)) = 0;

    temp = im2bw(temp,thresh);

    focal_adhesion_sections = temp | focal_adhesion_sections;

end

end