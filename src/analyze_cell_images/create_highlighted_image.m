function highlighted_image = create_highlighted_image(I,highlights,varargin)

image_size = size(I);

highlight_color_number = 2;
if (size(varargin,1) > 0)
    highlight_color_number = varargin{1};
end

if (size(image_size) < 3) 
    highlighted_image = zeros(image_size(1),image_size(2),3);
    highlighted_image(:,:,1) = I;
    highlighted_image(:,:,2) = I;
    highlighted_image(:,:,3) = I;
else
    highlighted_image = I;
end
    
temp_highlighted = I;

temp_highlighted(find(highlights)) = 1;

highlighted_image(:,:,highlight_color_number) = highlight_image;

end