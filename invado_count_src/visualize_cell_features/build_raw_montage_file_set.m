function build_raw_montage_file_set(base_dir,target_file)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('base_dir',@(x)exist(x,'dir') == 7);
i_p.addRequired('target_file',@(x)ischar(x));

i_p.parse(base_dir,target_file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

field = dir(base_dir);
%toss out first two entries they are '.' and '..', do same in below
%searches
field = field(3:end);

cell_composite = cell(0);
gel_composite = cell(0);

total_images = length(dir(fullfile(base_dir,field(1).name,'individual_pictures'))) - 2;

for i = 1:total_images    
    cell_images = cell(0);
    gel_images = cell(0);
    
    for j=1:length(field)
        cell_base = fullfile(base_dir,field(j).name,'Images/puncta');
        cell_files = dir(cell_base);
        cell_files = cell_files(3:end);
        cell_images{j} = double(imread(fullfile(cell_base,cell_files(i).name)));
        
        gel_base = fullfile(base_dir,field(j).name,'Images/gel');
        gel_files = dir(gel_base);
        gel_files = gel_files(3:end);
        gel_images{j} = double(imread(fullfile(gel_base,gel_files(i).name)));
    end
    
    %reverse every other column of images, to deal with the zig zag pattern    
    cell_temp = cell_images;
    gel_temp = gel_images;
    for j=1:(length(cell_images)/5)
        if (mod(j,2) == 1) 
            continue; 
        end
        for k=1:5
            cell_temp{j*5+k} = cell_images{j*5+(6-k)};
            gel_temp{j*5+k} = gel_images{j*5+(6-k)};
        end
    end
    cell_images = cell_temp;
    gel_images = gel_temp;
    
    cell_composite{i} = cat(1,cell_images{1:5});
    gel_composite{i} = cat(1,gel_images{1:5});
    
    for j=2:(length(cell_images)/5)
        cell_composite{i} = cat(2,cell_composite{i}, cat(1,cell_images{((j-1)*5+1):(j*5)}));
        gel_composite{i} = cat(2,gel_composite{i}, cat(1,gel_images{((j-1)*5+1):(j*5)}));
    end

    cell_composite{i} = imresize(cell_composite{i},[800,NaN]);
    gel_composite{i} = imresize(gel_composite{i},[800,NaN]);
    
    if (mod(i,10) == 0)
        disp(['Done processing: ',num2str(i)])
    end
end

cell_min_max = [min(min([cell_composite{1:end}])),max(max([cell_composite{1:end}]))];
gel_min_max = [min(min([gel_composite{1:end}])),max(max([gel_composite{1:end}]))];

for i=1:length(cell_composite)
    cell_composite{i} = (cell_composite{i} - cell_min_max(1))/range(cell_min_max);
    gel_composite{i} = (gel_composite{i} - gel_min_max(1))/range(gel_min_max);
end

if (not(exist(fullfile(output_dir,'cells'),'dir')))
    mkdir(fullfile(output_dir,'cells'));
end
if (not(exist(fullfile(output_dir,'gel'),'dir')))
    mkdir(fullfile(output_dir,'gel'));
end
if (not(exist(fullfile(output_dir,'both'),'dir')))
    mkdir(fullfile(output_dir,'both'));
end


for i=1:length(cell_composite)
    padding_length = num2str(length(num2str(total_images)));
    padded_i_num = sprintf(['%0',padding_length,'d'],i);
    
    imwrite(cell_composite{i}, fullfile(output_dir,'cells',[padded_i_num,'.png']));
    imwrite(gel_composite{i}, fullfile(output_dir,'gel',[padded_i_num,'.png']));
    imwrite(cat(3,gel_composite{i},cell_composite{i},zeros(size(cell_composite{i}))), fullfile(output_dir,'both',[padded_i_num,'.png']));

    if (mod(i,10) == 0)
        disp(['Done writing: ',num2str(i)])
    end    
end