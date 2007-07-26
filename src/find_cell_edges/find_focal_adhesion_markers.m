function focal_markers = find_focal_adhesion_markers(focal_image);

se = strel('disk',8);
Ie = imerode(focal_image, se);
Iobr = imreconstruct(Ie, focal_image);

Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);

focal_markers = imregionalmax(Iobrcbr);

end