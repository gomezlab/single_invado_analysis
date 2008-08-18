function high_image = create_highlighted_image(I,high,varargin)
%CREATE_HIGHLIGHTED_IMAGE    add highlights to an image
%
%   Color Number Codes:
%       -1: red
%       -2: green
%       -3: blue
%
%   H_I = create_highlighted_image(I,HIGHLIGHTS) adds green highlights to
%   image 'I', using the binary image 'HIGHLIGHTS' as the guide
%
%   H_I = create_highlighted_image(I,HIGHLIGHTS,'color',N) adds highlights of
%   color 'N' to image 'I', using the binary image 'HIGHLIGHTS' as the
%   guide
%
%   H_I = create_highlighted_image(I,HIGHLIGHTS,'color',[R,G,B]) adds
%   highlights of color specified by the RGB sequence '[R,G,B]' to image
%   'I', using the binary image 'HIGHLIGHTS' as the guide

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;
i_p.FunctionName = 'CREATE_HIGHLIGHTED_IMAGE';

i_p.addRequired('I',@(x)isnumeric(x) || islogical(x));
i_p.addRequired('high',@(x)(isnumeric(x) || islogical(x)));

i_p.parse(I,high);

i_p.addParamValue('color_map',[0,1,0],@(x)(all(high(:) == 0) || (isnumeric(x) && (size(x,1) == size(unique(high),1) - 1))));
i_p.addParamValue('mix_percent',1,@(x)(isnumeric(x)));

i_p.parse(I,high,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_size = size(I);

if (size(image_size) < 3)
    high_image = cat(3,I,I,I);
else
    high_image = I;
end

for j = 1:size(I,1)
    for k = 1:size(I,2)
        if (i_p.Results.high(j,k) == 0), continue; end

        this_cmap = i_p.Results.color_map(i_p.Results.high(j,k),:);
        
        if (all(this_cmap == 0)), continue; end
        
        assert(length(this_cmap) == size(high_image,3),'Error: wrong number of entries in color map');
        for l = 1:size(high_image,3)
            high_image(j,k,l) = this_cmap(l)*i_p.Results.mix_percent + high_image(j,k,l)*(1-i_p.Results.mix_percent);
        end
    end
end

end