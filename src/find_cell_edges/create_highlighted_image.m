function highlighted_image = create_highlighted_image(I,highlights,varargin)

image_size = size(I);

highlighted_image = zeros(image_size(1),image_size(2),3);
highlighted_image(:,:,1) = I;
highlighted_image(:,:,2) = I;
highlighted_image(:,:,3) = I;

highlight_color_number = 2;
if (size(varargin,1) > 0)
    highlight_color_number = varargin{1};
end

for i = 1:image_size(1)
    for j = 1:image_size(2)
        if (highlights(i,j))
            highlighted_image(i,j,highlight_color_number) = 1;
        end
    end
end

end