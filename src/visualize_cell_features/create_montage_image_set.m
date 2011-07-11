function montage=create_montage_image_set(image_set,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('image_set',@iscell);

i_p.addOptional('num_cols',0,@isnumeric);
i_p.addOptional('num_rows',0,@isnumeric);
i_p.addOptional('num_images',0,@isnumeric);
i_p.addOptional('phase','none',@ischar);
i_p.addOptional('pixel_size',-Inf,@(x)isnumeric(x) & x > 0);
i_p.addOptional('bar_size',10,@(x)isnumeric(x) & x > 0);
i_p.addOptional('bar_position',1,@(x)isnumeric(x) & x > 0);

i_p.parse(image_set,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The image set may start when many empty cells, first we clear all those
%out
while isempty(image_set{1}), image_set = image_set(2:end); end
while isempty(image_set{end}), image_set = image_set(1:(end-1)); end

if (isempty(strmatch('num_images', i_p.UsingDefaults)))
    assert(size(image_set,2) >= i_p.Results.num_images, ...
        'Problem with number of images requested: image_set:%d Requested Number:%d',size(image_set,2), i_p.Results.num_images)
    if (isempty(strmatch('phase', i_p.UsingDefaults)))
        if (strcmpi(i_p.Results.phase, 'assembly'))
            image_set = image_set(1:i_p.Results.num_images);
        elseif (strcmpi(i_p.Results.phase, 'disassembly'))
            image_set = image_set((size(image_set,2) - i_p.Results.num_images + 1):end);
        else
            warning('FA:phaseType','Expected assembly or disassembly for phase parameter, got %s',i_p.Results.phase)
        end
        assert(size(image_set,2) == i_p.Results.num_images, ...
            'Problem with removing images: image_set:%d Requested Number:%d',size(image_set,2), i_p.Results.num_images)
    else
        warning('FA:adCount','When num_images parameter specified, expected phase parameter to be set to either assembly or disassembly.')
    end
end

%convert all images in set to 3 layers, if not already three layers
for i = 1:size(image_set,2)
    if (length(size(image_set{i})) < 3)
        image_set{i} = cat(3,image_set{i},image_set{i},image_set{i});
    end
end

total_images = size(image_set,2);

%Determine how many images to put on each side
if (not(any(strcmp('num_cols',i_p.UsingDefaults))))
    images_per_side = [ceil(total_images/i_p.Results.num_cols),i_p.Results.num_cols];
elseif (not(any(strcmp('num_rows',i_p.UsingDefaults))))
    images_per_side = [i_p.Results.num_rows,ceil(total_images/i_p.Results.num_rows)];
else
    images_per_side = [floor(sqrt(total_images)),0];
    images_per_side(2) = ceil(total_images/images_per_side(1));
end

image_size = size(image_set{1});

montage = 0.5*ones(image_size(1)*images_per_side(1)+images_per_side(1) - 1, ...
    image_size(2)*images_per_side(2)+images_per_side(2) - 1, ...
    3);
for j = 1:images_per_side(1)
    for k = 1:images_per_side(2)
        i_index = (j-1)*images_per_side(2) + k;
        
        if (i_index > total_images), continue; end
        if (isempty(image_set{i_index})), continue; end
        
        montage((j-1)*image_size(1) + j:(j)*image_size(1) + j - 1, ...
            (k-1)*image_size(2) + k:(k)*image_size(2) + k - 1, ...
            1:3) = image_set{i_index};
    end
end

if (isempty(strmatch('pixel_size', i_p.UsingDefaults)))
    montage = draw_scale_bar(montage,i_p.Results.pixel_size, 'bar_size', ...
        i_p.Results.bar_size,'position_code',i_p.Results.bar_position);
end
