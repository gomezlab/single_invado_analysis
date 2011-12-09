base_dir = '../../data/hand_segmentations/control/';

folders = dir(base_dir);
folders = folders(3:end);
for i = 1:length(folders)
    folders(i).name = fullfile(base_dir,folders(i).name);
end

% for i = 1:length(folders)
%     determine_error_rates(folders(i).name);
%     disp(['Done with ', folders(i).name]);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consensus Stat Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(fullfile(folders(1).name,'consensus_stats.mat'));
area = consensus_props.Area;
 
for i = 2:length(folders)
    load(fullfile(folders(i).name,'consensus_stats.mat'));
    area = [area,consensus_props.Area];
end


% std_ranges = 0.5:0.5:12;
% min_size_vals = 1:10;
% 
% [x_vals, y_vals] = meshgrid(min_size_vals,std_ranges);
% 
% load(fullfile(folders(1).name,'stats.mat'));
% matches = stats.matches;
% comp_hits = stats.comp_hits;
% hand_hits = stats.hand_hits;
%  
% for i = 2:length(folders)
%     load(fullfile(folders(i).name,'stats.mat'));
%     matches = matches + stats.matches;
%     comp_hits = comp_hits + stats.comp_hits;
%     hand_hits = hand_hits + stats.hand_hits;
% end
% 
% false_neg = 100*(1 - matches./hand_hits);
% false_pos = 100*(1 - matches./comp_hits);
% 
% [c,h]=contourf(x_vals,y_vals,false_pos);
% xlabel('Puncta Seed Minimum (Pixels)','Fontsize',24);
% ylabel('Stdev Threshold','Fontsize',24);
% colorbar('FontSize',15);
% print('False_positives.pdf');
% clabel(c,h,'FontSize',15);
% 
% 
% [c,h]=contourf(x_vals,y_vals,false_neg);
% xlabel('Puncta Seed Minimum (Pixels)','Fontsize',24);
% ylabel('Stdev Threshold','Fontsize',24);
% colorbar('FontSize',15);
% print('False_negatives.pdf');
% clabel(c,h);

% plot(std_ranges,1-matches./hand_hits); ylim([0,1]); hold on; 
% plot(std_ranges,1 - matches./comp_hits,'r'); hold off;