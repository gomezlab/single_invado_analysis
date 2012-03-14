rm(list = ls())
source('../../../invado_count_src/analyze_cell_features/invado_analysis_lib.R')
source('misc_plotting_properties.R')
debug=TRUE;

################################################################################
# Data Reading/Initial Processing
################################################################################
base_dir = '../../../results/Invado_count/'

percents_files = list();
percents_files$Control = Sys.glob(file.path(base_dir,'*/*ontrol/*/cell_props/lin_time_series/Cell_gel_diff_percent.csv'))
percents_files$BB94 = Sys.glob(file.path(base_dir,'*/BB94/*/cell_props/lin_time_series/Cell_gel_diff_percent.csv'))
percents_files$DMSO = Sys.glob(file.path(base_dir,'*/DMSO/*/cell_props/lin_time_series/Cell_gel_diff_percent.csv'))
percents_files$FAK = Sys.glob(file.path(base_dir,'*/FAK/*/cell_props/lin_time_series/Cell_gel_diff_percent.csv'))
percents_files$PP2 = Sys.glob(file.path(base_dir,'*/PP2/*/cell_props/lin_time_series/Cell_gel_diff_percent.csv'))
percents_files$Noc = Sys.glob(file.path(base_dir,'*/Noc/*/cell_props/lin_time_series/Cell_gel_diff_percent.csv'))
percents_files$PurvA = Sys.glob(file.path(base_dir,'*/PurvA/*/cell_props/lin_time_series/Cell_gel_diff_percent.csv'))

percent_cumprods = list()
degrade_seq = list()

for (exp_type in names(percents_files)) {
    percent_cumprods[[exp_type]] = process_gel_diff_percents(percents_files[[exp_type]],min.lifetime=20)
    degrade_seq[[exp_type]] = percent_cumprods[[exp_type]]$seq

    percent_cumprods[[exp_type]]$seq <- NULL

    percent_cumprods[[exp_type]] = as.data.frame(percent_cumprods[[exp_type]])
    print(paste('Done reading', exp_type, 'got', length(percent_cumprods[[exp_type]]$cell_num)))
}


colors = get_exp_colors();

# source('../../../invado_count_src/analyze_cell_features/invado_analysis_lib.R')
# these_row_sums = get_row_cumsums(percents_files$Control,min.lifetime=20)

stop();
################################################################################
# Plotting
################################################################################

# temp = do.call(rbind, degrade_seq$BB94$degrade_seq_pad);
# 
# control_sub = subset(percent_cumprods$Control, min_degrade > -0.05)
# BB94_sub = subset(percent_cumprods$BB94, min_degrade > -0.05)
# BB94_sub = subset(percent_cumprods$BB94, max_degrade > 0.2)

# percent_changes = as.matrix(read.csv(percents_files$Control[3],header=F));
# percent_changes = percent_changes/100 + 1;
# 
# time = seq(1.5,by=0.5,along.with=percent_changes[1,])
# 
# svg('sample_cell_degrade.svg',width=3,height=3)
# par(bty='n',mar=c(2.6,2.5,0,0), mgp=c(1.6,0.5,0),xpd=T)
# plot(time[2:51], 1-cumprod(percent_changes[2,2:51]),typ='l',
#     xlab='Time Since Plating (hours)',ylab='Cummulative % Degradation',
#     ylim=c(0,0.7),xlim=c(0,max(time)),lwd=3,col='green');
# lines(time[2:51],1-cumprod(percent_changes[1,2:51]),typ='l',lwd=3,col='purple');
# # lines(time[3:51],cumprod(percent_changes[3,3:51]),typ='l',lwd=3,col='green');
# lines(time[3:51],1-cumprod(percent_changes[4,3:51]),typ='l',lwd=3);
# graphics.off()
