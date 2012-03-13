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

degrade_percent_distributions = list()

for (exp_type in names(percents_files)) {
    degrade_percent_distributions[[exp_type]] = get_non_nan_data(percents_files[[exp_type]],min.lifetime=20)
    print(paste('Done reading', exp_type, 'got', dim(degrade_percent_distributions[[exp_type]])[1]))
}

daily_percents = list()
hits = grep("2011_12_02", degrade_percent_distributions$Control$file)
daily_percents$control$day_1 = degrade_percent_distributions$Control[hits,]
hits = grep("2011_12_05", degrade_percent_distributions$Control$file)
daily_percents$control$day_2 = degrade_percent_distributions$Control[hits,]
hits = grep("2011_12_16", degrade_percent_distributions$Control$file)
daily_percents$control$day_3 = degrade_percent_distributions$Control[hits,]

colors = get_exp_colors();

percent_overlap_files = list();
percent_overlap_files$Control = Sys.glob(file.path(base_dir,'*/control/*/cell_props/lin_time_series/Overlap_percent.csv'))

overlap_percents = get_non_nan_data(percent_overlap_files$Control,min.lifetime=20)

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
    degrade_percent_distributions$DMSO$values,
    degrade_percent_distributions$FAK$values,
    degrade_percent_distributions$PP2$values,
    degrade_percent_distributions$Noc$values,
    degrade_percent_distributions$PurvA$values,
    outline=T,
    ylab='% Change in Fluorescence',axes=F)

axis(2)
axis(1,at=1:length(label_names),labels = label_names,padj=0.6)
graphics.off()

unique_BB94_misses = unique(subset(degrade_percent_distributions$BB94, values < -1.8)$file)

subset(degrade_percent_distributions$BB94, file == unique_BB94_misses[2] & values < -1.8)

svg('overlap_percents.svg',height=4,width=4)
par(bty='n',mar=c(2.7,2.5,0.45,0), mgp=c(1.6,0.5,0),xpd=T)
hist(overlap_percents$values,main='')
graphics.off()
