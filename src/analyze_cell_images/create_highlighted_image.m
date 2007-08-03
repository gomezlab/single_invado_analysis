function highlighted_image = create_highlighted_image(varargin)

%%Default Variable Values
highlight_color_number = 2;

%%Error Checking and Setup
if (size(varargin,2) < 2)
    error('ERROR: create_highlighted_image - not enough parameters specified, at least two needed');
elseif (isstruct(varargin{1}))
    image_data = varargin{1};
    if (isfield(image_data,'focal_edge_highlights'))
        I = image_data.focal_edge_highlights;
    else
        if (isfield(image_data,'focal_image'))
            I = image_data.focal_image;
        else
            error('ERROR: create_highlighted_image - Image data not found in struct');
        end
    end

    if (isnumeric(varargin{2}) || islogical(varargin{2}))
        highlights = varargin{2};
    elseif (isfield(image_data,varargin{2}))
        highlights = getfield(image_data,varargin{2});
    else
        error('ERROR: create_highlighted_image - Can''t find highlight matrix');
    end

    if (size(varargin,2) >= 3)
        if (isnumeric(varargin{3}) && varargin{3} < 4 && varargin{3} > 0)
            highlight_color_number = varargin{3};
        else
            error('ERROR: create_highlighted_image - 3rd parameter must be a number, indicating which color to highlight, 1-red, 2-green or 3-blue');
        end
    end
elseif (isnumeric(varargin{1}))
    I = varargin{1};

    if (isnumeric(varargin{2}) || islogical(varargin{2}))
        highlights = varargin{2};
    else
        error('ERROR: create_highlighted_image - If 1st parameter is numeric, so should the second');
    end

    if (size(varargin,2) >= 3)
        if (isnumeric(varargin{3}) && varargin{3} < 4 && varargin{3} > 0)
            highlight_color_number = varargin{3};
        else
            error('ERROR: create_highlighted_image - 3rd parameter must be a number, indicating which color to highlight, 1-red, 2-green or 3-blue');
        end
    end
else
    error('ERROR: create_highlighted_image - 1st parameter must be a struct or numeric array');
end

%%Main Program
image_size = size(I);

if (size(image_size) < 3)
    highlighted_image = zeros(image_size(1),image_size(2),3);
    highlighted_image(:,:,1) = I;
    highlighted_image(:,:,2) = I;
    highlighted_image(:,:,3) = I;
else
    highlighted_image = I;
end

if (size(image_size) < 3)
    temp_highlighted = I;
else
    temp_highlighted = I(:,:,highlight_color_number);
end

temp_highlighted(find(highlights)) = 1;

highlighted_image(:,:,highlight_color_number) = temp_highlighted;

end