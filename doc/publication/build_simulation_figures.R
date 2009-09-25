rm(list = ls())
source('FA_analysis_lib.R')
library(lattice)
library(geneplotter)
library(Hmisc)

################################################################################
#Result loading
################################################################################
raw_data <- list()
single_props <- list()

#Wild-type FA
exp_dirs <- Sys.glob('../../results/simulation/*/*/models')
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]

for (dir in exp_dirs) {
    this_data_type = basename(dirname(dirname(dir)))
    raw_data[[this_data_type]] = load_results(dir,file.path('intensity.Rdata'))
}

print('Done Loading Data')

########################################
#Result filtering
########################################
processed = list();

for (exp_type in names(raw_data)) {
    processed$no_filt[[exp_type]] = filter_results(raw_data[[exp_type]], 
        min_R_sq = -Inf, max_p_val = Inf);
    
    processed$only_signif[[exp_type]] = filter_results(raw_data[[exp_type]], 
        min_R_sq = -Inf, max_p_val = 0.05);
    
    processed$high_Rsq[[exp_type]] = filter_results(raw_data[[exp_type]]);
}

#rm(raw_data)
gc()

print('Done Filtering Data')
out_folder = '../../doc/publication/figures'
dir.create(out_folder,recursive=TRUE, showWarnings=FALSE);

################################################################################
#Plotting
################################################################################

########################################
#Stationary
########################################

longevity_filt = raw_data$stationary[[1]]$exp_props$longevity >= 25;
ad_sig_filt = raw_data$stationary[[1]]$exp_props$ad_sig[which(longevity_filt)]
mean_area_filt = raw_data$stationary[[1]]$exp_props$mean_area[which(longevity_filt)]

start_x_filt = raw_data$stationary[[1]]$exp_props$start_x[which(longevity_filt)]
x_clust = kmeans(start_x_filt, seq(min(start_x_filt), max(start_x_filt), length=10))
stopifnot(all(x_clust$size == 11 || x_clust$size == 10))

start_y_filt = raw_data$stationary[[1]]$exp_props$start_y[which(longevity_filt)]
y_clust = kmeans(start_y_filt, seq(min(start_y_filt), max(start_y_filt), length=11))
stopifnot(all(y_clust$size == 10 || y_clust$size == 9))

filt_data = data.frame(ad_sig = ad_sig_filt, mean_area = mean_area_filt, x_clusters = x_clust$cluster,
        y_clusters = y_clust$cluster)

dir.create(dirname(file.path(out_folder, 'simulation', 'stationary_hist.svg')), 
    recursive=TRUE, showWarnings=FALSE);
svg(file.path(out_folder, 'simulation', 'stationary_results.svg'), height=7*1.1, width=7*1.1);
par(mar=c(4,4,0.75,0), bty='n');
layout(rbind(c(1,2),c(3,4)))

plot.new()
mtext('A',adj=-.2,side=3,line=-0.75,cex=1.5);

hist_data = hist(raw_data$stationary[[1]]$exp_props$longevity, main='', xlab='Adhesion Longevity')
text(25, 300+hist_data$counts[which(hist_data$mids == 25)], sum(raw_data$stationary[[1]]$exp_props$longevity >= 25, na.rm=TRUE));
mtext('B',adj=-.2,side=3,line=-0.75,cex=1.5);

par(mar=c(4,4,1.5,0), bty='n');
boxplot(ad_sig ~ y_clusters, data=filt_data, notch=T, ylab='Average Adhesion Intensity', xlab='Column Number',
        names = 5:15)
mtext('C',adj=-.2,side=3,line=-0.75,cex=1.5);

boxplot(mean_area ~ x_clusters, data=filt_data, ylab = 'Mean Adhesion Area', xlab='Row Number')
mtext('D',adj=-.2,side=3,line=-0.75,cex=1.5);

graphics.off()

########################################
#Moving
########################################



########################################
#Phases
########################################

svg(file.path(out_folder, 'simulation', 'kinetics_results.svg'), height=7, width=7);
par(mar=c(4,4,0.75,0), bty='n');
layout(rbind(c(1,2),c(3,4)))

boxplot(processed$no_filt$phases_10$assembly$slope,
        processed$no_filt$phases_11$assembly$slope, 
        processed$no_filt$phases_12$assembly$slope, 
        processed$no_filt$phases_13$assembly$slope, 
        processed$no_filt$phases_14$assembly$slope, 
        processed$no_filt$phases_15$assembly$slope, 
        processed$no_filt$phases_16$assembly$slope, 
        processed$no_filt$phases_17$assembly$slope, 
        processed$no_filt$phases_18$assembly$slope, 
        processed$no_filt$phases_19$assembly$slope, 
        processed$no_filt$phases_20$assembly$slope,
        ylab = "Assembly Phase Slope", xlab = "Assembly Phase Lengths", names=10:20);
for (i in 1:11) {
        segments(i-0.5,log(0.4723663/0.2341950)/(i+9), i+0.5, log(0.4723663/0.2341950)/(i+9), col='red')
}
mtext('A',adj=-.2,side=3,line=-0.75,cex=1.5);

boxplot(processed$no_filt$phases_10$assembly$length,
        processed$no_filt$phases_11$assembly$length, 
        processed$no_filt$phases_12$assembly$length, 
        processed$no_filt$phases_13$assembly$length, 
        processed$no_filt$phases_14$assembly$length, 
        processed$no_filt$phases_15$assembly$length, 
        processed$no_filt$phases_16$assembly$length, 
        processed$no_filt$phases_17$assembly$length, 
        processed$no_filt$phases_18$assembly$length, 
        processed$no_filt$phases_19$assembly$length, 
        processed$no_filt$phases_20$assembly$length,
        ylab = "Detected Assembly Phase Length", xlab = "Assembly Phase Lengths", names=10:20);
# for (i in 1:11) {
#         segments(i-0.5,i+9, i+0.5, i+9, col='red')
# }
mtext('B',adj=-.2,side=3,line=-0.75,cex=1.5);

boxplot(processed$no_filt$phases_10$disassembly$slope,
        processed$no_filt$phases_11$disassembly$slope, 
        processed$no_filt$phases_12$disassembly$slope, 
        processed$no_filt$phases_13$disassembly$slope, 
        processed$no_filt$phases_14$disassembly$slope, 
        processed$no_filt$phases_15$disassembly$slope, 
        processed$no_filt$phases_16$disassembly$slope, 
        processed$no_filt$phases_17$disassembly$slope, 
        processed$no_filt$phases_18$disassembly$slope, 
        processed$no_filt$phases_19$disassembly$slope, 
        processed$no_filt$phases_20$disassembly$slope,
        ylab = "Assembly Phase Slope", xlab = "Assembly Phase Lengths", names=10:20);
for (i in 1:11) {
        segments(i-0.5,log(0.4723663/0.2341950)/(i+9), i+0.5, log(0.4723663/0.2341950)/(i+9), col='red')
}
mtext('C',adj=-.2,side=3,line=-0.75,cex=1.5);

boxplot(processed$no_filt$phases_10$disassembly$length,
        processed$no_filt$phases_11$disassembly$length, 
        processed$no_filt$phases_12$disassembly$length, 
        processed$no_filt$phases_13$disassembly$length, 
        processed$no_filt$phases_14$disassembly$length, 
        processed$no_filt$phases_15$disassembly$length, 
        processed$no_filt$phases_16$disassembly$length, 
        processed$no_filt$phases_17$disassembly$length, 
        processed$no_filt$phases_18$disassembly$length, 
        processed$no_filt$phases_19$disassembly$length, 
        processed$no_filt$phases_20$disassembly$length,
        ylab = "Detected Assembly Phase Length", xlab = "Assembly Phase Lengths", names=10:20);
# for (i in 1:11) {
#         segments(i-0.5,i+9, i+0.5, i+9, col='red')
# }
mtext('D',adj=-.2,side=3,line=-0.75,cex=1.5);
graphics.off()
