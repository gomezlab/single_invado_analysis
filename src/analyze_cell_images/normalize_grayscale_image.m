function [return_image] = normalize_grayscale_images(I)
% NORMALIZE_GRAYSCALE_IMAGE    using a given image, the same image with the
%                              minimum and maximum pixel value adjusted to
%                              0 and 1 respectively is returned in double
%                              format
%
%   normalize_grayscale_image(I) the image 'I' converted to a double and
%   the min and max values adjusted to 0 and 1 respectively


I = double(I);
image_size = size(I);

if (size(image_size,2) > 2) 
    if (image_size(3) > 1)
        I = I(:,:,1);
    end
end

I = I - min(I(:));
I = I/max(I(:));

return_image = I;

end