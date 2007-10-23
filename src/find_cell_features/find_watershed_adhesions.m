function focal_adhesion_sections = find_watershed_adhesions(varargin)
% FIND_WATERSHED_ADHESIONS    using a set of precomputed watershed
%                             boundaries, the adhesions within each
%                             boundary are extracted using a simple binary
%                             threshold
%
%   find_watershed_adhesions(I_struct) uses the focal adhesion protein
%   image 'I_struct.original_image' and the watershed edges identified in
%   'I_struct.watershed_edges' to locate focal adhesions, additionally, the
%   cell mask image must all be available in 'I_struct.cell_mask'
%
%   find_watershed_adhesions(I,wshed_edges,cell_mask) uses the focal
%   adhesion protein image 'I' and the watershed edges identified in
%   'wshed_edges' to locate focal adhesions, additionally, the cell mask
%   image must all be available in 'cell_mask'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (isempty(varargin))
    error('ERROR: find_watershed_adhesions - No parameters specified, see ''help find_watershed_adhesions''');
end

if (isstruct(varargin{1}))
    image_data = varargin{1};
    if (not(isfield(image_data,'original_image')))
        error('ERROR: find_watershed_adhesions - The first parameter is a struct, but the original focal image is not available in ''original_image''');
    end
    if (not(isfield(image_data,'cell_mask')))
        error('ERROR: find_watershed_adhesions - The first parameter is a struct, but the cell mask is not available in ''cell_mask''');
    end
    if (not(isfield(image_data,'watershed_edges')))
        error('ERROR: find_watershed_adhesions - The first parameter is a struct, but the watershed edges image is not available in ''watershed_edges''');
    end
    original_image = image_data.original_image;
    cell_mask = image_data.cell_mask;
    watershed_edges = image_data.watershed_edges;
else
    if (length(varargin) ~= 3)
        error('ERROR: find_watershed_adhesions - Not the correct number of parameters specified');
    else
        if (isnumeric(varargin{1}) && (isnumeric(varargin{2}) || islogical(varargin{2})) && (isnumeric(varargin{2}) || islogical(varargin{2})))
            original_image = varargin{1};
            watershed_edges = varargin{2};
            cell_mask = varargin{3};
        else
            error('ERROR: find_watershed_adhesions - Problem with one of the three parameters specified, check');
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
focal_adhesion_sections = zeros(size(original_image,1),size(original_image,2));
min_thresh = adaptive_thresh(original_image(find(cell_mask)));
watershed_labels = bwlabel(~watershed_edges.*cell_mask);

for i = 1:max(watershed_labels(:))
    thresh = adaptive_thresh(original_image(find(watershed_labels == i)));

    if (thresh < min_thresh)
        thresh = min_thresh;
    end
    temp = original_image;
    temp(find(watershed_labels ~= i)) = 0;

    temp = im2bw(temp,thresh);

    focal_adhesion_sections = temp | focal_adhesion_sections;

end

stats = regionprops(bwlabel(focal_adhesion_sections),'Area');
focal_adhesion_sections = ismember(bwlabel(focal_adhesion_sections),find([stats.Area] > 16));

end