function focal_markers = find_focal_adhesion_markers(varargin);
% FIND_FOCAL_ADHESION_MARKERS    locate the markers which will be used in
%                                the watershed segmentation
%
%   find_focal_adhesion_markers(image_data) locates the markers in
%   image 'image_data.original_image', using the binary image in
%   'image_data.cell_mask' to limit the possible locations for markers
%
%   find_focal_adhesion_markers(focal_image,cell_mask) locates the markers
%   in image 'focal_image', using the binary image in 'cell_mask' to limit
%   the possible locations for markers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Parameter and Variable Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (isempty(varargin))
    error('ERROR: find_focal_adhesion_markers - no parameters specified, see ''help find_focal_adhesion_markers''');
end

if (isstruct(varargin{1}))
    image_data = varargin{1};
    if (not(isfield(image_data,'original_image')))
        error('ERROR: find_focal_adhesion_markers - Missing focal image data in struct, looked in ''original_image''');
    end
    if (not(isfield(image_data,'cell_mask')))
        error('ERROR: find_focal_adhesion_markers - Missing cell mask data in struct, looked in ''cell_mask''');
    end

    focal_image = image_data.focal_image;
    cell_mask = image_data.cell_mask;
    
elseif (length(varargin) >= 2)
    if (not(isnumeric(varargin{1})))
        error('ERROR: find_focal_adhesion_markers - first parameter is not numeric');        
    elseif (not(isnumeric(varargin{2}) || islogical(varargin{2})))
        error('ERROR: find_focal_adhesion_markers - second parameter is not numeric or logical');
    end
    
    focal_image = varargin{1};
    cell_mask = varargin{2};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Command to sharpen the input image, note that unsharp actually sharpens
%the edges of the image

%focal_image = imfilter(focal_image,fspecial('unsharp'),'replicate');

thresh = adaptive_thresh(focal_image(find(cell_mask)));

focal_markers = im2bw(focal_image,thresh);
focal_markers = imfill(focal_markers,'holes');
focal_markers = bwmorph(focal_markers,'clean');

end