function [watershed_edges,varargout] = locate_watershed_edges(varargin)
% LOCATE_WATERSHED_EDGES    locates the watershed edges in a provided focal
%                           adhesion image
%
%   locate_watershed_edges(I_struct) using the focal image in
%   'I_struct.original_image' and the focal adhesion markers in
%   'I_struct.focal_markers', the watershed edges are forced to remain
%   bounded by the cell mask in 'I_stuct.cell_mask'
%
%   locate_watershed_edges(orig_I,f_m,c_m) using the focal image in
%   'orig_I' and the focal adhesion markers in 'f_m', the watershed edges
%   are forced to remain bounded by the cell mask in 'c_m'
%
%   NOTE: the binary image of the watershed edges are always returned from
%   this function, if an additional return value is requested, the labeled
%   watershed regions are returned along with the watershed edges


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (isempty(varargin))
    error('ERROR: locate_watershed_edges - Expected parameters, see ''help locate_watershed_edges''');
elseif (isstruct(varargin{1}))
    image_data = varargin{1};
    if (not(isfield(image_data,'original_image')))
        error('ERROR: locate_watershed_edges - with struct based parameters, expected to find original focal image in ''original_image''');
    elseif (not(isfield(image_data,'cell_mask')))
        error('ERROR: locate_watershed_edges - with struct based parameters, expected to find cell mask image in ''cell_mask''');
    elseif (not(isfield(image_data,'focal_markers')))
        error('ERROR: locate_watershed_edges - with struct based parameters, expected to find focal markers image in ''focal_markers''');
    end

    original_image = image_data.original_image;
    focal_markers = image_data.focal_markers;
    cell_mask = image_data.cell_mask;
else
    if (length(varargin) == 3)
        if (not(isnumeric(varargin{1})))
            error('ERROR: locate_watershed_edges - The first parameter should be the original image, a numeric');
        elseif (not(isnumeric(varargin{2}) || islogical(varargin{2})))
            error('ERROR: locate_watershed_edges - The second parameter should be the watershed_edges, a numeric or logical');
        elseif (not(isnumeric(varargin{3}) || islogical(varargin{3})))
            error('ERROR: locate_watershed_edges - The third parameter should be the cell_edge, a numeric or logical');
        end
    else
        error('ERROR: locate_watershed_edges - When not using struct parameters, expected three parameters');
    end

    original_image = varargin{1};
    focal_markers = varargin{2};
    cell_mask = varargin{3};
end

%the cell edge will be needed following the watershed segmentation
cell_edge = bwperim(cell_mask);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inverted_focal_image = -original_image + 1;
inverted_focal_image = imfilter(inverted_focal_image,fspecial('gaussian'));
inverted_focal_image(~cell_mask) = 1.1;
inverted_focal_image = imimposemin(inverted_focal_image,focal_markers);

watershed_labels = watershed(inverted_focal_image);
watershed_labels(~cell_mask) = 0;
watershed_labels(cell_edge) = 0;

watershed_edges = zeros(size(original_image,1),size(original_image,2));
watershed_edges(find(watershed_labels == 0)) = 1;
watershed_edges(~cell_mask) = 0;

if (nargout > 0)
    varargout{1} = watershed_labels;
end