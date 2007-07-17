function [return_image] = normalize_grayscale_images(I)

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