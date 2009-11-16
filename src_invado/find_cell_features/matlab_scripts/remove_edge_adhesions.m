function cleaned_binary = remove_edge_adhesions(threshed_image)

edge_border = zeros(size(threshed_image));
edge_border(1,:) = 1;
edge_border(end,:) = 1;
edge_border(:,1) = 1;
edge_border(:,end) = 1;

[row_indexes,col_indexes] = ind2sub(size(threshed_image), find(edge_border));
edge_adhesions = bwselect(threshed_image,col_indexes,row_indexes,4);

cleaned_binary = threshed_image & not(edge_adhesions);
