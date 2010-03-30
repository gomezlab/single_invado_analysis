rm(list = ls())
source('../../src_invado/analyze_cell_features/invado_analysis_lib.R')
source('../../src/analyze_cell_features/FA_analysis_lib.R')
library(lattice)
library(geneplotter)
library(Hmisc)
library(boot)

debug=TRUE;

################################################################################
# Data Reading/Initial Processing
################################################################################

dirs_to_load = list();
dirs_to_load$control = Sys.glob('../../results/Invadopodia/control/*/adhesion_props/');
dirs_to_load$coro1C_kd = Sys.glob('../../results/Invadopodia/Coro1C/*/adhesion_props/');
dirs_to_load$coro1B_kd = Sys.glob('../../results/Invadopodia/Coro1B/*/adhesion_props/');
dirs_to_load$CotL_kd = Sys.glob('../../results/Invadopodia/CotL/*/adhesion_props/');

raw_data = list();
for (exp_type in names(dirs_to_load)) {
    raw_data[[exp_type]] = load_results_data_frame(dirs_to_load[[exp_type]], 
        file.path('models','puncta_props_corr.Rdata'), 'all_props');
}

data_sets = list();
for (i in names(raw_data)) {
    invado_name = paste("invado_", i, sep='');
    non_invado_name = paste("non_invado_", i, sep='');
    
    filter_sets = build_filter_sets(raw_data[[i]]);
    
    data_sets[[invado_name]] = subset(raw_data[[i]], filter_sets$invado_filter);
    
    data_sets[[non_invado_name]] = subset(raw_data[[i]], filter_sets$non_invado_filter);
}

data_sets$day_1_control = subset(data_sets$invado_control, regexpr('11_16',data_sets$invado_control$experiment) == 1)
data_sets$day_2_control = subset(data_sets$invado_control, regexpr('11_22',data_sets$invado_control$experiment) == 1)
data_sets$day_3_control = subset(data_sets$invado_control, regexpr('02_23',data_sets$invado_control$experiment) == 1)

data_sets$day_1_coro1C = subset(data_sets$invado_coro1C_kd, regexpr('11_14',data_sets$invado_coro1C_kd$experiment) == 1)
data_sets$day_2_coro1C = subset(data_sets$invado_coro1C_kd, regexpr('11_15',data_sets$invado_coro1C_kd$experiment) == 1)
data_sets$day_3_coro1C = subset(data_sets$invado_coro1C_kd, regexpr('11_22',data_sets$invado_coro1C_kd$experiment) == 1)

ts_props = list()

for (exp_type in names(dirs_to_load)) {
    ts_props[[exp_type]]$bleaching_curve = load_data_files(dirs_to_load[[exp_type]], c('bleaching_curves.csv'), debug=F, inc_exp_names=T);
}

stop()

################################################################################
# Plotting
################################################################################

########################################
# Exploritory
########################################

####################
# Area Comparisons
####################

area_sets_to_use = list(data_sets$invado_control$largest_area, 
                        data_sets$invado_coro1B_kd$largest_area, 
                        data_sets$invado_coro1C_kd$largest_area, 
                        data_sets$invado_CotL_kd$largest_area)

data_names = c('Control','Coro1B KD', 'Coro1C KD', 'Coactosin KD')

svg('area_boxplots.svg')
boxplot(area_sets_to_use, names=data_names,
        ylab='Largest Puncta Area (\u03BCm\u00B2)')
graphics.off()

area_conf_ints = gather_barplot_properties(area_sets_to_use);

svg('area_barplots.svg')
par(mar=c(2,4,0.5,0))
x_pos = barplot(area_conf_ints$mean, names=data_names,
        ylab='Largest Puncta Area (\u03BCm\u00B2)', ylim=c(0,max(area_conf_ints$yplus)))
errbar(t(x_pos),area_conf_ints$mean,area_conf_ints$yplus, area_conf_ints$yminus,
       add=TRUE,cex=1E-10,lwd=1.5)
graphics.off()

####################
# Longevity Comparisons
####################

longev_sets_to_use = list(data_sets$invado_control$longevity*5, 
                        data_sets$invado_coro1B_kd$longevity*5, 
                        data_sets$invado_coro1C_kd$longevity*5, 
                        data_sets$invado_CotL_kd$longevity*5)

data_names = c('Control','Coro1B KD', 'Coro1C KD', 'Coactosin KD')

svg('longevity_boxplots.svg')
boxplot(longev_sets_to_use, names=data_names,
        ylab='Longevity (min)')
graphics.off()

longev_conf_ints = gather_barplot_properties(longev_sets_to_use);

svg('longevity_barplots.svg')
par(mar=c(2,4,0.5,0))
x_pos = barplot(longev_conf_ints$mean, names=data_names,
        ylab='Longevity (min)', ylim=c(0,max(longev_conf_ints$yplus)))
errbar(t(x_pos),longev_conf_ints$mean,longev_conf_ints$yplus, longev_conf_ints$yminus,
       add=TRUE,cex=1E-10,lwd=1.5)
graphics.off()

########################################
# Day to Day Comparisons
########################################

####################
# Longevity Comparisons
####################

data_sets_to_use = list(data_sets$day_1_control$longevity*5, 
                        data_sets$day_2_control$longevity*5, 
                        data_sets$day_3_control$longevity*5);

data_names = c('Day 1','Day 2', 'Day 3')

boxplot(data_sets_to_use, names=data_names,
        ylab='Longevity (min)')

longev_conf_ints = gather_barplot_properties(data_sets_to_use);

par(mar=c(2,4,0.5,0))
x_pos = barplot(longev_conf_ints$mean, names=data_names,
        ylab='Longevity (min)', ylim=c(0,max(longev_conf_ints$yplus)))
errbar(t(x_pos),longev_conf_ints$mean,longev_conf_ints$yplus, longev_conf_ints$yminus,
       add=TRUE,cex=1E-10,lwd=1.5)


boxplot(data_sets$day_1_control$mean_vals, data_sets$day_2_control$mean_vals,data_sets$day_3_control$mean_vals)
boxplot(data_sets$day_1_control_uncorr$mean_vals, data_sets$day_2_control_uncorr$mean_vals,data_sets$day_3_control_uncorr$mean_vals)

boxplot_names = c('Ctrl Day 1', 'Coro1C Day 1', 'Coro1C Day 1', 'Ctrl Day 2', 'Coro1C Day 2', 'Ctrl Day 3', 'CotL Day 3')
boxplot(data_sets$day_1_control$mean_vals, data_sets$day_1_coro1C$mean_vals, data_sets$day_2_coro1C$mean_vals, 
    data_sets$day_2_control$mean_vals, data_sets$day_3_coro1C$mean_vals, 
    data_sets$day_3_control$mean_vals, data_sets$invado_CotL_kd$mean_vals, notch=T, names=boxplot_names)

boxplot(data_sets$day_1_control$largest_area, data_sets$day_1_coro1C$largest_area, 
    data_sets$day_2_coro1C$largest_area, data_sets$day_2_control$largest_area, 
    data_sets$day_3_coro1C$largest_area, data_sets$day_3_control$largest_area, 
    data_sets$invado_CotL_kd$largest_area, 
    notch=T, names=boxplot_names)

boxplot(data_sets$day_1_control$longevity, data_sets$day_1_coro1C$longevity, data_sets$day_2_coro1C$longevity, 
    data_sets$day_2_control$longevity, data_sets$day_3_coro1C$longevity, 
    data_sets$day_3_control$longevity, data_sets$invado_CotL_kd$longevity, notch=T, names=boxplot_names)

bleach_curves_mat = matrix(NA,ncol = length(ts_props$bleaching_curve[[1]]), nrow = 121);
percent_diffs = c();
for (i in 1:length(ts_props$bleaching_curve[[1]])) {
    d_s_length = length(ts_props$bleaching_curve[[1]][[i]][1,]);
    bleach_curves_mat[1:d_s_length,i] = as.numeric(ts_props$bleaching_curve[[1]][[i]][1,]);

    this_percent = min(ts_props$bleaching_curve[[1]][[i]][1,])/max(ts_props$bleaching_curve[[1]][[i]][1,]);
    percent_diffs = c(percent_diffs, this_percent);
}

correlations = list();
for (i in 1:length(ts_props$bleaching_curve[[1]])) {
    correlations$overall_vs_notcell = c(correlations$overall_vs_notcell,
        cor(as.numeric(ts_props$bleaching_curve[[1]][[i]][1,]),as.numeric(ts_props$bleaching_curve[[1]][[i]][2,])));
    
    correlations$overall_vs_puncta = c(correlations$overall_vs_puncta,
        cor(as.numeric(ts_props$bleaching_curve[[1]][[i]][1,]),as.numeric(ts_props$bleaching_curve[[1]][[i]][3,])));

    correlations$notcell_vs_puncta = c(correlations$notcell_vs_puncta,
        cor(as.numeric(ts_props$bleaching_curve[[1]][[i]][2,]),as.numeric(ts_props$bleaching_curve[[1]][[i]][3,])));
}

props_to_look_at = c('longevity','mean_vals', 'largest_area');
layout(cbind(c(1,2),c(3,4)))
for (i in props_to_look_at) {
    min_val = min(c(data_sets$invado_control[[i]], data_sets$invado_coro1B_kd[[i]], data_sets$invado_coro1C_kd[[i]]));
    if (i == 'mean_vals') {
        boxplot(data_sets$invado_control_corr[[i]], data_sets$invado_coro1B_kd_corr[[i]], data_sets$invado_coro1C_kd_corr[[i]], data_sets$invado_CotL_kd_corr[[i]], notch=T, main=i, names=c('ctrl','1B','1C', 'CotL'), ylim=c(min_val,0))
    } else {
        boxplot(data_sets$invado_control_corr[[i]], data_sets$invado_coro1B_kd_corr[[i]], data_sets$invado_coro1C_kd_corr[[i]], data_sets$invado_CotL_kd_corr[[i]], notch=T, main=i, names=c('ctrl','1B','1C', 'CotL'))
    }
}

layout(cbind(1,2,3))

########################################
# Longevity Bi-hist
########################################
bin_size = 50;
max_data = 5*max(data_sets$invado_control$longevity, data_sets$invado_coro1C_kd$longevity);
max_data = ceil(max_data/bin_size)*bin_size;
min_data = 5*min(data_sets$invado_control$longevity, data_sets$invado_coro1C_kd$longevity);

data_1 = hist(data_sets$invado_control$longevity*5, breaks=seq(min_data,max_data,by=bin_size));
data_2 = hist(data_sets$invado_coro1C_kd$longevity*5, breaks=data_1$breaks);
graphics.off()

all_data = rbind(data_1$counts/sum(data_1$counts), data_2$counts/sum(data_2$counts));
range_names = c()
for (i in 1:(length(data_1$breaks)-1)) {
    range_names = c(range_names, paste(data_1$breaks[i],"-",data_1$breaks[i+1],sep=''));
}

svg('longevity_hist.svg')
par(mar = c(5.75,4,0.25,0))
barplot(all_data, beside=T, names = range_names, ylab="Percentage", 
        #las controls the orientation of the labels, las=2 is perpendicular to axis
        legend = c('Control', 'Coronin1C KD'), las=2)
mtext("Longevity Ranges (min)", side=1, line=4.5)
graphics.off()

plot(data_sets$invado_control$longevity*5, data_sets$invado_control$largest_area, pch=20, cex=0.75)
points(data_sets$invado_coro1C_kd$longevity*5, data_sets$invado_coro1C_kd$largest_area, col='red', pch=20, cex=0.75)

t.test(data_sets$invado_control$longevity*5, data_sets$invado_coro1C_kd$longevity*5)
t.test(data_sets$invado_control$largest_area, data_sets$invado_coro1C_kd$largest_area)

########################################
# Largest Size Bi-hist
########################################
bin_size = 0.5;
max_data = max(data_sets$invado_control$largest_area, data_sets$invado_coro1C_kd$largest_area);
max_data = ceil(max_data/bin_size)*bin_size;
min_data = min(data_sets$invado_control$largest_area, data_sets$invado_coro1C_kd$largest_area);
min_data = floor(min_data/bin_size)*bin_size;

data_1 = hist(data_sets$invado_control$largest_area, breaks=seq(min_data,max_data,by=bin_size));
data_2 = hist(data_sets$invado_coro1C_kd$largest_area, breaks=data_1$breaks);
graphics.off()

all_data = rbind(data_1$counts/sum(data_1$counts), data_2$counts/sum(data_2$counts));
range_names = c()
for (i in 1:(length(data_1$breaks)-1)) {
    range_names = c(range_names, paste(data_1$breaks[i],"-",data_1$breaks[i+1],sep=''));
}

svg('largest_area_hist.svg')
par(mar = c(4.75,4,0.25,0))
barplot(all_data, beside=T, names = range_names, ylab="Percentage", 
        #las controls the orientation of the labels, las=2 is perpendicular to axis
        legend = c('Control', 'Coronin1C KD'), las=2)
mtext('Max Invadopodia Area (\u03BCm\u00B2)', side=1, line=3.5)
graphics.off()

########################################
# Average Local Diff Values Bi-hist
########################################
bin_size = 0.5E-3*10^2;

control_mean_vals = data_sets$invado_control_corr$mean_vals*-10^2;
coro1C_kd_mean_vals = data_sets$invado_coro1C_kd$mean_vals*-10^2;
coro1B_kd_mean_vals = data_sets$invado_coro1B_kd_corr$mean_vals*-10^2;

max_data = max(control_mean_vals, coro1C_kd_mean_vals);
max_data = ceil(max_data/bin_size)*bin_size;
min_data = min(control_mean_vals, coro1C_kd_mean_vals);
min_data = floor(min_data/bin_size)*bin_size;

data_1 = hist(control_mean_vals, breaks=seq(min_data,max_data,by=bin_size));
data_2 = hist(coro1C_kd_mean_vals, breaks=data_1$breaks);
data_3 = hist(coro1B_kd_mean_vals, breaks=data_1$breaks);
graphics.off()

all_data = rbind(data_1$counts/sum(data_1$counts), data_2$counts/sum(data_2$counts))
#, data_3$counts/sum(data_3$counts));
range_names = c()
for (i in 1:(length(data_1$breaks)-1)) {
    range_names = c(range_names, paste(data_1$breaks[i],"-",data_1$breaks[i+1],sep=''));
}

svg('mean_vals_hist.svg')
par(mar = c(6,4,0.25,0))
barplot(all_data, beside=T, names = range_names, ylab="Percentage", 
        #las controls the orientation of the labels, las=2 is perpendicular to axis
        las=2)
mtext('Mean Local Difference (arbitrary units * -10\u00B2)', side=1, line=4.75)
legend('topright', c('Control', 'Coronin1C KD'), fill = c('gray30','gray90'))
graphics.off()

plot(data_sets$invado_control$longevity, data_sets$invado_control$mean_vals, ylim=c(min_data,max_data))
points(data_sets$invado_coro1C_kd$longevity, data_sets$invado_coro1C_kd$mean_vals, col='red')

