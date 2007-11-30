function [watershed_edges,varargout] = locate_watershed_edges(original_image,focal_markers,cell_mask)
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

i_p = inputParser;
i_p.FunctionName = 'LOCATE_WATERSHED_EDGES';

i_p.addRequired('original_image',@isnumeric);
i_p.addRequired('focal_markers',@(x) isnumeric(x) || islogical(x));
i_p.addRequired('cell_mask',@(x) isnumeric(x) || islogical(x));

i_p.parse(original_image,focal_markers,cell_mask);

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