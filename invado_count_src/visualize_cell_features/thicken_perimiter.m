function perim_thick = thicken_perimiter(perims,objects,varargin)

if (islogical(perims)), perims = double(perims); end
if (islogical(objects)), objects = double(objects); end

thickness = 3;
if (size(varargin,2) > 0)
    thickness = varargin{1};
end

perim_thick = perims;
three_disk = strel('disk',thickness,0);
for cell_num = 1:max(perim_thick(:))
    this_cell_perim = perims == cell_num;
    this_cell = objects == cell_num;
    
    this_cell_perim = imdilate(this_cell_perim,three_disk);
    this_cell_perim = this_cell_perim & this_cell;
    
    perim_thick(this_cell_perim) = cell_num;
end
