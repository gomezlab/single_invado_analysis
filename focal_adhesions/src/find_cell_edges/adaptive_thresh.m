function [threshold] = adpative_thresh(I)

threshold = 0.5;

threshold_change_thresh = 0.0001;

upper_mean = mean(I(find(I > threshold)));
lower_mean = mean(I(find(I <= threshold)));

new_thresh = 0.2*upper_mean + 0.8*lower_mean;

while (abs(threshold - new_thresh) > threshold_change_thresh)
    threshold = new_thresh;
    
    upper_mean = mean(I(find(I > threshold)));
    lower_mean = mean(I(find(I <= threshold)));

    new_thresh = 0.2*upper_mean + 0.8*lower_mean;
end

