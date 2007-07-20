function highlighted_image = create_highlighted_image(I,highlights,varargin)

image_size = size(I);

highlight_color_number = 2;
if (size(varargin,1) > 0)
    highlight_color_number = varargin{1};
end

highlighted_image = zeros(image_size(1),image_size(2),3);
highlighted_image(:,:,1) = I;
highlighted_image(:,:,2) = I;
highlighted_image(:,:,3) = I;

highlight_pixels = find(highlights);

highlight_image = I;
for i = 1:size(highlight_pixels,1)
    highlight_image(highlight_pixels(i)) = 1;
end

highlighted_image(:,:,highlight_color_number) = highlight_image;

end