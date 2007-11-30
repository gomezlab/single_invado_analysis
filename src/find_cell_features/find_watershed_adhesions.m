function focal_adhesion_sections = find_watershed_adhesions(I,watershed_edges,cell_mask)
% FIND_WATERSHED_ADHESIONS    using a set of precomputed watershed
%                             boundaries, the adhesions within each
%                             boundary are extracted using a simple binary
%                             threshold
%
%   find_watershed_adhesions(I,wshed_edges,cell_mask) uses the focal
%   adhesion protein image 'I' and the watershed edges identified in
%   'wshed_edges' to locate focal adhesions, additionally, the cell mask
%   image must all be available in 'cell_mask'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;
i_p.FunctionName = 'FIND_WATERSHED_ADHESIONS';

i_p.addRequired('I',@isnumeric);
i_p.addRequired('watershed_edges',@(x)isnumeric(x) || islogical(x));
i_p.addRequired('cell_mask',@(x)isnumeric(x) || islogical(x));

i_p.parse(I,watershed_edges,cell_mask);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
focal_adhesion_sections = zeros(size(I,1),size(I,2));
min_thresh = adaptive_thresh(I(cell_mask));
watershed_labels = bwlabel(~watershed_edges.*cell_mask);

for i = 1:max(watershed_labels(:))
    thresh = adaptive_thresh(I(watershed_labels == i));

    if (thresh < min_thresh)
        thresh = min_thresh;
    end
    temp = I;
    temp(find(watershed_labels ~= i)) = 0;

    temp = im2bw(temp,thresh);

    focal_adhesion_sections = temp | focal_adhesion_sections;

end

end