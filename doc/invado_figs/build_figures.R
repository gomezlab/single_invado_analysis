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

dir.create('area',recursive=TRUE, showWarnings=FALSE);

area_sets_to_use = list(data_sets$invado_control$mean_area, 
                        data_sets$invado_coro1B_kd$mean_area, 
                        data_sets$invado_coro1C_kd$mean_area, 
                        data_sets$invado_CotL_kd$mean_area)

data_names = c('Control','Coro1B KD', 'Coro1C KD', 'Coactosin KD')

for (i in 1:length(data_names)) {
    data_names[[i]] = paste(data_names[[i]], " (n=", length(na.omit(area_sets_to_use[[i]])), ")", sep='');
}

svg(file.path('area','area_boxplots.svg'), width=9)
boxplot(area_sets_to_use, names=data_names,
        ylab='Mean Puncta Area (\u03BCm\u00B2)')
graphics.off()

area_conf_ints = gather_barplot_properties(area_sets_to_use);

svg(file.path('area','area_barplots.svg'), width=8)
par(mar=c(2,4,0.5,0))
x_pos = barplot(area_conf_ints$mean, names=data_names,
        ylab='Mean Puncta Area (\u03BCm\u00B2)', ylim=c(0,max(area_conf_ints$yplus)))
errbar(t(x_pos),area_conf_ints$mean,area_conf_ints$yplus, area_conf_ints$yminus,
       add=TRUE,cex=1E-10,lwd=1.5)
graphics.off()

####################
# Longevity Comparisons
####################

dir.create('longevity',recursive=TRUE, showWarnings=FALSE);

longev_sets_to_use = list(data_sets$invado_control$longevity*5, 
                          data_sets$invado_coro1B_kd$longevity*5, 
                          data_sets$invado_coro1C_kd$longevity*5, 
                          data_sets$invado_CotL_kd$longevity*5)

svg(file.path('longevity','longevity_boxplots.svg'), width=9)
boxplot(longev_sets_to_use, names=data_names,
        ylab='Longevity (min)')
graphics.off()

longev_conf_ints = gather_barplot_properties(longev_sets_to_use);

svg(file.path('longevity','longevity_barplots.svg'), width=8)
par(mar=c(2,4,0.5,0))
x_pos = barplot(longev_conf_ints$mean, names=data_names,
        ylab='Longevity (min)', ylim=c(0,max(longev_conf_ints$yplus)))
errbar(t(x_pos),longev_conf_ints$mean,longev_conf_ints$yplus, longev_conf_ints$yminus,
       add=TRUE,cex=1E-10,lwd=1.5)
graphics.off()

####################
# Longevity/Area Comparisons
####################

dir.create('longevity_size_comp',recursive=TRUE, showWarnings=FALSE);

data_names = c('Control','Coro1B KD', 'Coro1C KD', 'Coactosin KD')

for (i in 1:length(longev_sets_to_use)) {
    svg(file.path('longevity_size_comp',paste(data_names[[i]], '.svg', sep='')))
    par(bty='n');
    plot(longev_sets_to_use[[i]], area_sets_to_use[[i]], main=data_names[[i]],
        xlab='Longevity (min)', ylab='Mean Puncta Area (\u03BCm\u00B2)')
    graphics.off()
}

svg(file.path('longevity_size_comp','all.svg'))
par(bty='n');
layout(rbind(c(1,2),c(3,4)))
for (i in 1:length(longev_sets_to_use)) {
    plot(longev_sets_to_use[[i]], area_sets_to_use[[i]], main=data_names[[i]],
        xlab='Longevity (min)', ylab='Mean Puncta Area (\u03BCm\u00B2)')
}
graphics.off()

########################################
# Day to Day Comparisons
########################################

####################
# Area Comparisons
####################

dir.create('day_to_day',recursive=TRUE, showWarnings=FALSE);

area_sets_to_use = list(
        data_sets$day_1_control$mean_area,
        data_sets$day_2_control$mean_area,
        data_sets$day_3_control$mean_area,
        data_sets$invado_control$mean_area, 
        data_sets$invado_coro1B_kd$mean_area, 
        data_sets$invado_coro1C_kd$mean_area, 
        data_sets$invado_CotL_kd$mean_area)

data_names = c('Ctrl Day 1','Ctrl Day 2','Ctrl Day 3', 'Control','Coro1B KD', 'Coro1C KD', 'Coactosin KD')

for (i in 1:length(data_names)) {
    data_names[[i]] = paste(data_names[[i]], " (n=", length(na.omit(area_sets_to_use[[i]])), ")", sep='');
}

svg(file.path('day_to_day','area_boxplots.svg'), width=15)
boxplot(area_sets_to_use, names=data_names,
        ylab='Mean Puncta Area (\u03BCm\u00B2)')
graphics.off()

area_conf_ints = gather_barplot_properties(area_sets_to_use);

svg(file.path('day_to_day','area_barplots.svg'), width=15)
par(mar=c(2,4,0.5,0))
x_pos = barplot(area_conf_ints$mean, names=data_names,
        ylab='Mean Puncta Area (\u03BCm\u00B2)', ylim=c(0,max(area_conf_ints$yplus)))
errbar(t(x_pos),area_conf_ints$mean,area_conf_ints$yplus, area_conf_ints$yminus,
       add=TRUE,cex=1E-10,lwd=1.5)
graphics.off()

####################
# Longevity Comparisons
####################

dir.create('day_to_day',recursive=TRUE, showWarnings=FALSE);

longev_sets_to_use = list(
        data_sets$day_1_control$longevity*5,
        data_sets$day_2_control$longevity*5,
        data_sets$day_3_control$longevity*5,
        data_sets$invado_control$longevity*5, 
        data_sets$invado_coro1B_kd$longevity*5, 
        data_sets$invado_coro1C_kd$longevity*5, 
        data_sets$invado_CotL_kd$longevity*5
        )

svg(file.path('day_to_day','longevity_boxplots.svg'), width=15)
boxplot(longev_sets_to_use, names=data_names,
        ylab='Longevity (min)')
graphics.off()

longev_conf_ints = gather_barplot_properties(longev_sets_to_use);

svg(file.path('day_to_day','longevity_barplots.svg'), width=15)
par(mar=c(2,4,0.5,0))
x_pos = barplot(longev_conf_ints$mean, names=data_names,
        ylab='Longevity (min)', ylim=c(0,max(longev_conf_ints$yplus)))
errbar(t(x_pos),longev_conf_ints$mean,longev_conf_ints$yplus, longev_conf_ints$yminus,
       add=TRUE,cex=1E-10,lwd=1.5)
graphics.off()

########################################
# Bleaching Comparisons
########################################

bleach_curves_mat = matrix(NA, ncol = length(ts_props$coro1B_kd$bleaching_curve[[1]]), nrow = length(ts_props$coro1B_kd$bleaching_curve[[1]][[1]][,1]))

percent_diffs = c();
for (i in 1:length(ts_props$coro1B_kd$bleaching_curve[[1]])) {
    
    this_curve = ts_props$coro1B_kd$bleaching_curve[[1]][[i]][,1];

    bleach_curves_mat[,i] = this_curve;
    
    percent_diffs = c(percent_diffs, min(this_curve)/max(this_curve));
}


# bleach_curves_mat = matrix(NA,ncol = length(ts_props$bleaching_curve[[1]]), nrow = 121);
# percent_diffs = c();
# for (i in 1:length(ts_props$bleaching_curve[[1]])) {
#     d_s_length = length(ts_props$bleaching_curve[[1]][[i]][1,]);
#     bleach_curves_mat[1:d_s_length,i] = as.numeric(ts_props$bleaching_curve[[1]][[i]][1,]);
# 
#     this_percent = min(ts_props$bleaching_curve[[1]][[i]][1,])/max(ts_props$bleaching_curve[[1]][[i]][1,]);
#     percent_diffs = c(percent_diffs, this_percent);
# }

# correlations = list();
# for (i in 1:length(ts_props$bleaching_curve[[1]])) {
#     correlations$overall_vs_notcell = c(correlations$overall_vs_notcell,
#         cor(as.numeric(ts_props$bleaching_curve[[1]][[i]][1,]),as.numeric(ts_props$bleaching_curve[[1]][[i]][2,])));
#     
#     correlations$overall_vs_puncta = c(correlations$overall_vs_puncta,
#         cor(as.numeric(ts_props$bleaching_curve[[1]][[i]][1,]),as.numeric(ts_props$bleaching_curve[[1]][[i]][3,])));
# 
#     correlations$notcell_vs_puncta = c(correlations$notcell_vs_puncta,
#         cor(as.numeric(ts_props$bleaching_curve[[1]][[i]][2,]),as.numeric(ts_props$bleaching_curve[[1]][[i]][3,])));
# }
