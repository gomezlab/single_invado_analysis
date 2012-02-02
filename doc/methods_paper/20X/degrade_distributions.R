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

# degrade_percents$FAK = Sys.glob(file.path(base_dir,'*/FAK/overall*/degrade*'))
# degrade_percents$PP2 = Sys.glob(file.path(base_dir,'*/PP2/overall*/degrade*'))
# degrade_percents$Noc = Sys.glob(file.path(base_dir,'*/Noc/overall*/degrade*'))
# degrade_percents$PurvA = Sys.glob(file.path(base_dir,'*/PurvA/overall*/degrade*'))

degrade_percent_distributions = list()

for (exp_type in names(percents_files)) {
    degrade_percent_distributions[[exp_type]] = get_non_nan_data(percents_files[[exp_type]],min.lifetime=20)
    print(paste('Done reading', exp_type, 'got', dim(degrade_percent_distributions[[exp_type]])[1]))
}

mean_files = list();
mean_files$Control = Sys.glob(file.path(base_dir,'*/*ontrol/*/cell_props/lin_time_series/Cell_gel_diff_median.csv'))
mean_files$BB94 = Sys.glob(file.path(base_dir,'*/BB94/*/cell_props/lin_time_series/Cell_gel_diff_median.csv'))
mean_files$DMSO = Sys.glob(file.path(base_dir,'*/DMSO/*/cell_props/lin_time_series/Cell_gel_diff_median.csv'))

# degrade_mean$FAK = Sys.glob(file.path(base_dir,'*/FAK/overall*/degrade*'))
# degrade_mean$PP2 = Sys.glob(file.path(base_dir,'*/PP2/overall*/degrade*'))
# degrade_mean$Noc = Sys.glob(file.path(base_dir,'*/Noc/overall*/degrade*'))
# degrade_mean$PurvA = Sys.glob(file.path(base_dir,'*/PurvA/overall*/degrade*'))

degrade_mean_distributions = list()

for (exp_type in names(mean_files)) {
    degrade_mean_distributions[[exp_type]] = get_non_nan_data(mean_files[[exp_type]],min.lifetime=20)
    print(paste('Done reading', exp_type, 'got', dim(degrade_mean_distributions[[exp_type]])[1]))
}

colors = get_exp_colors();

stop();
################################################################################
# Plotting
################################################################################

label_names = c()
for (exp_type in names(degrade_percent_distributions)) {
    label_names = c(label_names,
        paste(exp_type,'\n(n=',length(degrade_percent_distributions[[exp_type]]$values),')',sep=''))
}

svg(file.path('degradation_distributions.svg'),width=3.5,height=3.5);
par(bty='n',mar=c(2.7,2.5,0.45,0), mgp=c(1.6,0.5,0),xpd=T)
boxplot(degrade_percent_distributions$Control$values,
    degrade_percent_distributions$BB94$values,
    degrade_percent_distributions$DMSO$values,outline=FALSE,
    ylab='% Change in Fluorescence',axes=F)

axis(2)
axis(1,at=1:3,labels = label_names,padj=0.6)
graphics.off()

unique_BB94_misses = unique(subset(degrade_percent_distributions$BB94, values < -0.85)$file)

subset(degrade_percent_distributions$BB94, file == unique_BB94_misses[3] & values < -0.85)
