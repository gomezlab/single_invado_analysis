function high_image = create_highlighted_image(I,high,varargin)
%CREATE_HIGHLIGHTED_IMAGE    add highlights to an image
%
%   Color Number Codes:
%       -1: red
%       -2: green
%       -3: blue
%
%   H_I = create_highlighted_image(I,HIGHLIGHTS) adds green highlights to
%   image 'I', using the binary image 'HIGHLIGHTS' as the guide
%
%   H_I = create_highlighted_image(I,HIGHLIGHTS,'color',N) adds highlights of 
%   color 'N' to image 'I', using the binary image 'HIGHLIGHTS' as the
%   guide
%
%   H_I = create_highlighted_image(I,HIGHLIGHTS,'color',[R,G,B]) adds
%   highlights of color specified by the RGB sequence '[R,G,B]' to image
%   'I', using the binary image 'HIGHLIGHTS' as the guide 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;
i_p.FunctionName = 'CREATE_HIGHLIGHTED_IMAGE';

i_p.addRequired('I',@(x)isnumeric(x) || islogical(x));
i_p.addRequired('high',@(x)(isnumeric(x) || islogical(x)));

i_p.addParamValue('color',2,@(x)(isnumeric(x) && (length(x) == 1 || length(x) == 3)));

i_p.parse(I,high,varargin{:});

if (length(i_p.Results.color) == 1)
    if (i_p.Results.color == 1)
        c_map = [1,0,0];
    elseif (i_p.Results.color == 2)
        c_map = [0,1,0];
    elseif (i_p.Results.color == 3)
        c_map = [0,0,1];
    end
else
    c_map = i_p.Results.color;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_size = size(I);

if (size(image_size) < 3)
    high_image = zeros(image_size(1),image_size(2),3);
    high_image(:,:,1) = I;
    high_image(:,:,2) = I;
    high_image(:,:,3) = I;
else
    high_image = I;
end

if(size(high) > 1)
    high = high(:,:,1);
end

pix_count = image_size(1)*image_size(2);

for i = 1:3
    high_image(find(high)+(i-1)*pix_count) = c_map(i);
end

end