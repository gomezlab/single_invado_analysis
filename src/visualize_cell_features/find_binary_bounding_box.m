function bbox = find_binary_bounding_box(I)
% FIND_BINARY_BOUNDING_BOX    finds the minimum box which completely encloses
%                             the provided binary image
%
%   [min_row,min_col,max_row,max_col] = find_binary_bounding_box(I) finds the min
%   and max, x and y coordinates that completely enclose the binary 
%   image, 'I'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.addRequired('I',@(x)isnumeric(x) || islogical(x));

i_p.parse(I);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp = regionprops(bwlabel(I),'BoundingBox');

bounds = zeros(size(temp,1),4);

for i = 1:size(temp,1)
    bounds(i,:) = temp(i).BoundingBox;
end

bbox = [Inf, Inf, -Inf, -Inf];

for i = 1:size(bounds,1)
    x1 = ceil(bounds(i,1));
    y1 = ceil(bounds(i,2));
    x2 = floor(x1 + bounds(i,3));
    y2 = floor(y1 + bounds(i,4));

    if (bbox(1) > x1), bbox(1) = x1; end
    if (bbox(2) > y1), bbox(2) = y1; end
    if (bbox(3) < x2), bbox(3) = x2; end
    if (bbox(4) < y2), bbox(4) = y2; end
end

if (bbox(1) < 0), bbox(1) = 0; end
if (bbox(2) < 0), bbox(2) = 0; end
if (bbox(3) > size(I,2)), bbox(3) = size(I,2); end
if (bbox(4) > size(I,2)), bbox(4) = size(I,2); end