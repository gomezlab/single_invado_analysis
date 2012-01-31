clear all;

files = dir('combined_results');
files = files(3:end);

all_series = cell(0);

for i = 1:length(files)
    load(fullfile('combined_results',files(i).name));
    fields = fieldnames(results);
    for j=1:length(fields)
        if (not(isfield(all_series,fields{j})))
            all_series.(fields{j}) = results.(fields{j});
        else
            all_series.(fields{j}) = [all_series.(fields{j}),results.(fields{j})];
        end
    end
end

[true_pos,false_neg,t,AUC_mean_no_corr] = perfcurve(all_series.degrade_or_not,all_series.diff_mean_no_corr,0);

plot(true_pos,false_neg);
hold on; 
plot([0,1],[0,1]);

[true_pos,false_neg,t,AUC_mean] = perfcurve(all_series.degrade_or_not,all_series.diff_mean,0);
plot(true_pos,false_neg,'r');

[true_pos,false_neg,t,AUC_median] = perfcurve(all_series.degrade_or_not,all_series.diff_median,0);
plot(true_pos,false_neg,'g');

[true_pos,false_neg,t,AUC_median] = perfcurve(all_series.degrade_or_not,all_series.diff_median,0);
plot(true_pos,false_neg,'g');

[true_pos,false_neg,t,AUC_percent] = perfcurve(all_series.degrade_or_not,all_series.diff_percent,0);
plot(true_pos,false_neg,'k');

hold off;

% [true_pos,false_neg,t,AUC_median] = perfcurve(all_series.degrade_or_not,all_series.diff_median,0);
% plot(true_pos,false_neg,'g');