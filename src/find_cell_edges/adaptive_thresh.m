function threshold = adpative_thresh(I,varargin)

upper_threshold_weight = 0.5;
if (size(varargin,1) > 0)
    upper_threshold_weight = varargin{1};
end

threshold = mean(I(:));

threshold_change_thresh = 0.0001;

upper_mean = mean(I(find(I > threshold)));
lower_mean = mean(I(find(I <= threshold)));

new_thresh = upper_threshold_weight*upper_mean + (1-upper_threshold_weight)*lower_mean;

while (abs(threshold - new_thresh) > threshold_change_thresh)
    threshold = new_thresh;
    
    upper_mean = mean(I(find(I > threshold)));
    lower_mean = mean(I(find(I <= threshold)));

    new_thresh = upper_threshold_weight*upper_mean + (1-upper_threshold_weight)*lower_mean;
end

end