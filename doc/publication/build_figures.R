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
exp_dirs <- Sys.glob('../../results/focal_adhesions/*/adhesion_props/models/')
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
raw_data$results = load_results(exp_dirs,file.path('intensity.Rdata'))
raw_data$corr_results = load_results(exp_dirs,file.path('local_corrected.Rdata'))
raw_data$CB_results = load_results(exp_dirs,file.path('CB_corrected.Rdata'))
raw_data$area = load_results(exp_dirs,file.path('area.Rdata'))
raw_data$ind_results <- load_data_files(exp_dirs, file.path('..','individual_adhesions.csv'), headers=T, debug=FALSE, inc_exp_names=FALSE);

background_correlation = load_results_2(exp_dirs, 'background_corr.Rdata')
single_props$fa$not_ad_int <- load_data_files(exp_dirs,
	file.path('..','single_props','Cell_not_ad_mean_intensity.csv'), 
	headers=F, debug=FALSE, inc_exp_names=FALSE);
single_props$fa$cell_int <- load_data_files(exp_dirs,
	file.path('..','single_props','Cell_mean_intensity.csv'), 
	headers=F, debug=FALSE, inc_exp_names=FALSE);
single_props$fa$ad_int <- load_data_files(exp_dirs,
	file.path('..','single_props','Adhesion_mean_intensity.csv'), 
	headers=F, debug=FALSE, inc_exp_names=FALSE);
single_props$fa$not_ad_int <- load_data_files(exp_dirs,
	file.path('..','single_props','Cell_not_ad_mean_intensity.csv'), 
	headers=F, debug=FALSE, inc_exp_names=FALSE);
single_props$fa$outside <- load_data_files(exp_dirs,
	file.path('..','single_props','Outside_mean_intensity.csv'), 
	headers=F, debug=FALSE, inc_exp_names=FALSE);

exp_dirs_reduced <- Sys.glob('../../results/focal_adhesions_reduced/*/adhesion_props/models/')
exp_dirs_reduced <- exp_dirs_reduced[file_test('-d',exp_dirs_reduced)]
raw_data$results_reduced = load_results(exp_dirs_reduced,file.path('intensity.Rdata'))

#S178A Results
exp_dirs_S <- Sys.glob('../../results/S178A/*/adhesion_props/models/')
exp_dirs_S <- exp_dirs_S[file_test('-d',exp_dirs_S)]
raw_data$results_S = load_results(exp_dirs_S,file.path('intensity.Rdata'))
raw_data$corr_results_S = load_results(exp_dirs_S,file.path('local_corrected.Rdata'))
raw_data$CB_results_S = load_results(exp_dirs_S,file.path('CB_corrected.Rdata'))
raw_data$area_S = load_results(exp_dirs_S,file.path('area.Rdata'))
raw_data$ind_results_S <- load_data_files(exp_dirs_S,file.path('..','individual_adhesions.csv'), headers=T, debug=FALSE, inc_exp_names=FALSE);

background_correlation_S = load_results_2(exp_dirs_S, 'background_corr.Rdata')
single_props$S178A = list()
single_props$S178A$cell_int <- load_data_files(exp_dirs_S,
	file.path('..','single_props','Cell_mean_intensity.csv'), 
	headers=F, debug=FALSE, inc_exp_names=FALSE);
single_props$S178A$not_ad_int <- load_data_files(exp_dirs_S,
	file.path('..','single_props','Cell_not_ad_mean_intensity.csv'), 
	headers=F, debug=FALSE, inc_exp_names=FALSE);
single_props$S178A$ad_int <- load_data_files(exp_dirs_S,
	file.path('..','single_props','Adhesion_mean_intensity.csv'), 
	headers=F, debug=FALSE, inc_exp_names=FALSE);
single_props$S178A$outside <- load_data_files(exp_dirs_S,
	file.path('..','single_props','Outside_mean_intensity.csv'), 
	headers=F, debug=FALSE, inc_exp_names=FALSE);

exp_dirs_S <- Sys.glob('../../results/S178A_reduced/*/adhesion_props/models/')
exp_dirs_S <- exp_dirs_S[file_test('-d',exp_dirs_S)]
raw_data$results_S_reduced = load_results(exp_dirs_S,file.path('intensity.Rdata'))

print('Done Loading Data')

########################################
#Result filtering
########################################
results_props = gather_general_props(raw_data$results)
results_S_props = gather_general_props(raw_data$results_S)

ind_exp_filt = gather_single_image_props(raw_data$ind_results)
ind_exp_filt_S = gather_single_image_props(raw_data$ind_results_S)

results_nofilt = filter_results(raw_data$results, min_R_sq = -Inf, max_p_val = Inf, cell_intensities=single_props$fa$cell_int)
results_reduced_nofilt = filter_results(raw_data$results_reduced, min_R_sq = -Inf, max_p_val = Inf)
results_S_nofilt = filter_results(raw_data$results_S, min_R_sq = -Inf, max_p_val = Inf, cell_intensities=single_props$S178A$cell_int)
results_S_reduced_nofilt = filter_results(raw_data$results_S_reduced, min_R_sq = -Inf, max_p_val = Inf, cell_intensities=single_props$S178A$cell_int)

corr_results_nofilt = filter_results(raw_data$corr_results, min_R_sq = -Inf, max_p_val = Inf, cell_intensities=single_props$fa$cell_int)

results_onlysignif = filter_results(raw_data$results, min_R_sq = -Inf, max_p_val = 0.05, cell_intensities=single_props$fa$cell_int)
results_S_onlysignif = filter_results(raw_data$results_S, min_R_sq = -Inf, max_p_val = 0.05, cell_intensities=single_props$S178A$cell_int)
corr_results_onlysignif = filter_results(raw_data$corr_results, min_R_sq = -Inf, max_p_val = 0.05)
corr_results_S_onlysignif = filter_results(raw_data$corr_results_S, min_R_sq = -Inf, max_p_val = 0.05)
CB_results_onlysignif = filter_results(raw_data$CB_results, min_R_sq = -Inf, max_p_val = 0.05)
CB_results_S_onlysignif = filter_results(raw_data$CB_results_S, min_R_sq = -Inf, max_p_val = 0.05)

area_onlysignif = filter_results(raw_data$area, min_R_sq = -Inf, max_p_val = 0.05)
area_filt = filter_results(raw_data$area)
area_S_onlysignif = filter_results(raw_data$area_S, min_R_sq = -Inf, max_p_val = 0.05)
results_filt = filter_results(raw_data$results, cell_intensities=single_props$fa$cell_int)
results_reduced_filt = filter_results(raw_data$results_reduced)
results_S_filt = filter_results(raw_data$results_S, cell_intensities=single_props$S178A$cell_int)
results_S_reduced_filt = filter_results(raw_data$results_S_reduced, cell_intensities=single_props$S178A$cell_int)

#rm(raw_data)
gc()

print('Done Filtering Data')
out_folder = '../../doc/publication/figures'
dir.create(out_folder,recursive=TRUE, showWarnings=FALSE);

stop()

assembly_nums = intersect(results_filt$assembly$lin_num, results_filt$disassembly$lin_num)

################################################################################
#Plotting
################################################################################
out_folder = '../../doc/publication/figures'
dir.create(out_folder,recursive=TRUE, showWarnings=FALSE);

########################################
#Control Checks
########################################

#Assembly and Disassembly Rates versus Expression
dir.create(dirname(file.path(out_folder,'controls','rates_vs_expression.svg')), 
    recursive=TRUE, showWarnings=FALSE);
svg(file.path(out_folder,'controls','rates_vs_expression.svg'))
layout(rbind(c(1,2), c(3,4)))

par(bty='n', mar=c(4,4.5,2,0))

exp_num_counts = c()
for (i in sort(unique(results_onlysignif$assembly$exp_num))) {
        exp_num_counts = c(exp_num_counts, length(which(results_onlysignif$assembly$exp_num == i)));
}

#Assembly/Disassembly Rates versus mean cell intensity
ylim_vals = range(c(results_onlysignif$control$assembly_slope, results_onlysignif$control$disassembly_slope))
plot(results_onlysignif$control$mean_cell_int, results_onlysignif$control$assembly_slope, 
	col='green', ylim=ylim_vals, xlab='Mean Cellular Paxillin Intensity', 
	ylab=expression(paste('Mean Rate (',min^-1,')',sep='')), main='Wild-type', pch=19);
points(results_onlysignif$control$mean_cell_int, results_onlysignif$control$disassembly_slope, col='red', pch=19)

legend('topright', c('Assembly','Disassembly'), fill=c('green','red'))
summary(lm(results_onlysignif$control$assembly_slope ~ results_onlysignif$control$mean_cell_int, weights=exp_num_counts))
summary(lm(results_onlysignif$control$disassembly_slope ~ results_onlysignif$control$mean_cell_int))
mtext('A',adj=-.2,side=3,line=-0.5,cex=1.5);

#birth/death position versus mean cell intensity
ylim_vals = range(c(results_onlysignif$control$birth_pos, results_onlysignif$control$death_pos), na.rm=T)
plot(results_onlysignif$control$mean_cell_int, results_onlysignif$control$birth_pos, 
	col='green', ylim=ylim_vals, xlab='Mean Cellular Paxillin Intensity', 
	ylab=expression(paste('Mean Distance from Edge (',mu,'m)',sep='')), main='Wild-type', pch=19);
points(results_onlysignif$control$mean_cell_int, results_onlysignif$control$death_pos, col='red', pch=19)

legend('topright', c('Birth','Death'), fill=c('green','red'))
summary(lm(results_onlysignif$control$birth_pos ~ results_onlysignif$control$mean_cell_int, weights=exp_num_counts))
summary(lm(results_onlysignif$control$death_pos ~ results_onlysignif$control$mean_cell_int, weights=exp_num_counts))
mtext('B',adj=-.2,side=3,line=-0.5,cex=1.5)

ylim_vals = range(c(results_S_onlysignif$control$assembly_slope, results_S_onlysignif$control$disassembly_slope), na.rm=T)
plot(results_S_onlysignif$control$mean_cell_int, results_S_onlysignif$control$assembly_slope, 
	col='green', ylim=ylim_vals, xlab='Mean Cellular Paxillin Intensity', 
	ylab=expression(paste('Mean Rate (',min^-1,')',sep='')), main='S178A', pch=19);
points(results_S_onlysignif$control$mean_cell_int, results_S_onlysignif$control$disassembly_slope, col='red', pch=19)

legend('topright', c('Assembly','Disassembly'), fill=c('green','red'))
summary(lm(results_S_onlysignif$control$assembly_slope ~ results_S_onlysignif$control$mean_cell_int))
summary(lm(results_S_onlysignif$control$disassembly_slope ~ results_S_onlysignif$control$mean_cell_int))
mtext('C',adj=-.2,side=3,line=-0.5,cex=1.5);

ylim_vals = range(c(results_S_onlysignif$control$birth_pos, results_S_onlysignif$control$death_pos), na.rm=T)
plot(results_S_onlysignif$control$mean_cell_int, results_S_onlysignif$control$birth_pos, 
	col='green', ylim=ylim_vals, xlab='Mean Cellular Paxillin Intensity', 
	ylab=expression(paste('Mean Distance from Edge (',mu,'m)',sep='')), main='S178A', pch=19);
points(results_S_onlysignif$control$mean_cell_int, results_S_onlysignif$control$death_pos, col='red', pch=19)


legend('topright', c('Birth','Death'), fill=c('green','red'))
summary(lm(results_S_onlysignif$control$birth_pos ~ results_S_onlysignif$control$mean_cell_int))
summary(lm(results_S_onlysignif$control$death_pos ~ results_S_onlysignif$control$mean_cell_int))
mtext('D',adj=-.2,side=3,line=-0.5,cex=1.5)
graphics.off()

#Rates versus expression unfiltered
svg(file.path(out_folder,'controls','rates_vs_expression_unfilt.svg'))
layout(rbind(c(1,2), c(3,4)))

par(bty='n', mar=c(4,4.5,2,0))

ylim_vals = range(c(results_onlysignif$control$assembly_slope, results_onlysignif$control$disassembly_slope))
plot(results_onlysignif$control$mean_cell_int, results_onlysignif$control$assembly_slope, 
	col='green', ylim=ylim_vals, xlab='Mean Paxillin Intensity', 
	ylab=expression(paste('Rate (',min^-1,')',sep='')), main='Wild-type', pch=19);
points(results_onlysignif$control$mean_cell_int, results_onlysignif$control$disassembly_slope, col='red', pch=19)
legend('topright', c('Assembly','Disassembly'), fill=c('green','red'))
summary(lm(results_onlysignif$control$assembly_slope ~ results_onlysignif$control$mean_cell_int))
summary(lm(results_onlysignif$control$disassembly_slope ~ results_onlysignif$control$mean_cell_int))
mtext('A',adj=-.2,side=3,line=-0.5,cex=1.5);

ylim_vals = range(c(results_onlysignif$control$birth_pos, results_onlysignif$control$death_pos), na.rm=T)
plot(results_onlysignif$control$mean_cell_int, results_onlysignif$control$birth_pos, 
	col='green', ylim=ylim_vals, xlab='Mean Paxillin Intensity', 
	ylab=expression(paste('Rate (',min^-1,')',sep='')), main='Wild-type', pch=19);
points(results_onlysignif$control$mean_cell_int, results_onlysignif$control$death_pos, col='red', pch=19)
legend('topright', c('Birth','Death'), fill=c('green','red'))
mtext('B',adj=-.2,side=3,line=-0.5,cex=1.5)

ylim_vals = range(c(results_S_filt$control$assembly_slope, results_S_filt$control$disassembly_slope), na.rm=T)
plot(results_S_filt$control$mean_cell_int, results_S_filt$control$assembly_slope, 
	col='green', ylim=ylim_vals, xlab='Mean Paxillin Intensity', 
	ylab=expression(paste('Rate (',min^-1,')',sep='')), main='S178A', pch=19);
points(results_S_filt$control$mean_cell_int, results_S_filt$control$disassembly_slope, col='red', pch=19)
legend('topright', c('Assembly','Disassembly'), fill=c('green','red'))
mtext('C',adj=-.2,side=3,line=-0.5,cex=1.5);

ylim_vals = range(c(results_S_filt$control$birth_pos, results_S_filt$control$death_pos), na.rm=T)
plot(results_S_filt$control$mean_cell_int, results_S_filt$control$birth_pos, 
	col='green', ylim=ylim_vals, xlab='Mean Paxillin Intensity', 
	ylab=expression(paste('Rate (',min^-1,')',sep='')), main='S178A', pch=19);
points(results_S_filt$control$mean_cell_int, results_S_filt$control$death_pos, col='red', pch=19)
legend('topright', c('Birth','Death'), fill=c('green','red'))
mtext('D',adj=-.2,side=3,line=-0.5,cex=1.5)
graphics.off()

#Resampled Slope
svg(file.path(out_folder,'controls','resampled_slopes.svg'))
layout(rbind(c(1,2),c(3,4)))
par(bty='n', mar=c(3,4.2,0,0));

boxplot_with_points(list(results_nofilt$assembly$slope, results_reduced_nofilt$assembly$slope/2), names=c('All', 'Sampled'), notch=T, ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')))
plot_size = par("usr");
text(0.9*plot_size[2], 0.9*plot_size[4], sprintf('%0.3f',t.test(results_nofilt$assembly$slope, results_reduced_nofilt$assembly$slope/2)$p.value), col='red')

boxplot_with_points(list(results_nofilt$disassembly$slope, results_reduced_nofilt$disassembly$slope/2), names=c('All', 'Sampled'), notch=T, ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')))
plot_size = par("usr");
text(0.9*plot_size[2], 0.9*plot_size[4], sprintf('%0.3f',t.test(results_nofilt$disassembly$slope, results_reduced_nofilt$disassembly$slope/2)$p.value), col='red')

boxplot_with_points(list(results_nofilt$assembly$R_sq, results_reduced_nofilt$assembly$R_sq), names=c('All', 'Sampled'), notch=T, ylab=expression(paste('Assembly ',R^2,sep='')))
plot_size = par("usr");
text(0.9*plot_size[2], 0.9*plot_size[3], sprintf('%0.3f',t.test(results_nofilt$assembly$R_sq, results_reduced_nofilt$assembly$R_sq)$p.value), col='red')

boxplot_with_points(list(results_nofilt$disassembly$R_sq, results_reduced_nofilt$disassembly$R_sq), names=c('All', 'Sampled'), notch=T, ylab=expression(paste('Disassembly ',R^2,sep='')))
plot_size = par("usr");
text(0.9*plot_size[2], 0.9*plot_size[3], sprintf('%0.3f',t.test(results_nofilt$disassembly$R_sq, results_reduced_nofilt$disassembly$R_sq)$p.value), col='red')
graphics.off()

#Reduced S178A
svg(file.path(out_folder,'controls','resampled_slopes_S178A.svg'))
layout(rbind(c(1,2),c(3,4)))
par(bty='n', mar=c(3,4.2,0,0));

boxplot_with_points(list(results_S_onlysignif$assembly$slope, results_S_reduced_nofilt$assembly$slope/2), names=c('All', 'Sampled'), notch=T, ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')))
plot_size = par("usr");
text(0.9*plot_size[2], 0.9*plot_size[4], sprintf('%0.3f',t.test(results_S_nofilt$assembly$slope, results_S_reduced_nofilt$assembly$slope/2)$p.value), col='red')

boxplot_with_points(list(results_S_nofilt$disassembly$slope, results_S_reduced_nofilt$disassembly$slope/2), names=c('All', 'Sampled'), notch=T, ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')))
plot_size = par("usr");
text(0.9*plot_size[2], 0.9*plot_size[4], sprintf('%0.3f',t.test(results_S_nofilt$disassembly$slope, results_S_reduced_nofilt$disassembly$slope/2)$p.value), col='red')

boxplot_with_points(list(results_S_nofilt$assembly$R_sq, results_S_reduced_nofilt$assembly$R_sq), names=c('All', 'Sampled'), notch=T, ylab=expression(paste('Assembly ',R^2,sep='')))
plot_size = par("usr");
text(0.9*plot_size[2], 0.9*plot_size[3], sprintf('%0.3f',t.test(results_S_nofilt$assembly$R_sq, results_S_reduced_nofilt$assembly$R_sq)$p.value), col='red')

boxplot_with_points(list(results_S_nofilt$disassembly$R_sq, results_S_reduced_nofilt$disassembly$R_sq), names=c('All', 'Sampled'), notch=T, ylab=expression(paste('Disassembly ',R^2,sep='')))
plot_size = par("usr");
text(0.9*plot_size[2], 0.9*plot_size[3], sprintf('%0.3f',t.test(results_S_nofilt$disassembly$R_sq, results_S_reduced_nofilt$disassembly$R_sq)$p.value), col='red')
graphics.off()

#Poor R^2 adhesion sequence plots
pdf(file.path(out_folder,'controls','poor_assembly.pdf'),width=14, height=7)
layout(cbind(c(1),c(2)))
for (i in which(results_nofilt$assembly$R_sq < 0.2)) {
    plot_ad_seq(raw_data$results[[results_nofilt$assembly$exp_num[i]]], results_nofilt$assembly$lin_num[i], main=results_nofilt$assembly$exp_dir[i], sub=results_nofilt$assembly$lin_num[i])
    plot_ad_seq(raw_data$area[[results_nofilt$assembly$exp_num[i]]], results_nofilt$assembly$lin_num[i], main=results_nofilt$assembly$exp_dir[i], sub=results_nofilt$assembly$lin_num[i])
}
graphics.off()

pdf(file.path(out_folder,'controls','poor_disassembly.pdf'),width=14, height=7)
layout(cbind(c(1),c(2)))
for (i in which(results_nofilt$disassembly$R_sq < 0.2)) {
    plot_ad_seq(raw_data$results[[results_nofilt$disassembly$exp_num[i]]], results_nofilt$disassembly$lin_num[i], main=results_nofilt$disassembly$exp_dir[i], sub=results_nofilt$disassembly$lin_num[i], type='disassembly')
    plot_ad_seq(raw_data$area[[results_nofilt$disassembly$exp_num[i]]], results_nofilt$disassembly$lin_num[i], main=results_nofilt$disassembly$exp_dir[i], sub=results_nofilt$disassembly$lin_num[i], type='disassembly')
}
graphics.off()

########################################
#All by all plots
########################################
#svg(file.path(out_folder,'as_splom.svg'),width=12,height=12)
#splom(cbind(results_filt$a[1:10]), sub='Assembly Scatter Plot (Filtered)')
#splom(cbind(results_nofilt$a[1:10]), sub='Assembly Scatter Plot')
#graphics.off()
#
#svg(file.path(out_folder,'dis_splom.svg'),width=12,height=12)
#splom(cbind(results_filt$dis[1:10]), sub='Disassembly Scatter Plot (Filtered)')
#splom(cbind(results_nofilt$dis[1:10]), sub='Disassembly Scatter Plot')
#graphics.off()


########################################
#Statics Properties
########################################
dir.create(dirname(file.path(out_folder,'statics','statics.svg')), 
    recursive=TRUE, showWarnings=FALSE)
svg(file.path(out_folder,'statics','statics.svg'),height=8)

layout_mat = rbind(c(rep(1,4),rep(2,4),rep(3,4)),
                   c(rep(4,6),rep(5,6)),
                   c(rep(6,6),rep(7,6)),
                   c(rep(0,3), rep(8,6), rep(0, 3)))
layout(layout_mat,heights=c(1,0.65,0.65,0.65))
par(bty='n', mar=c(0,4,1.6,0))

#Place holders for cell images
plot.new()
mtext('A',adj=-.31,side=3,line=0,cex=1.5)

plot.new()
mtext('B',adj=-.31,side=3,line=0,cex=1.5)

plot.new()
mtext('C',adj=-.31,side=3,line=0,cex=1.5)

#Histograms
par(bty='n', mar=c(4,4.2,0,0.1))
area_data = ind_exp_filt$Area[ind_exp_filt$Area < 5];
area_hist = hist(area_data, main="", ylab = "FA Count", 
	 xlab = expression(paste('FA Area (', symbol("m"), m^2, ')',sep='')));
    
#area_hist_model = lm(log(area_hist$counts) ~area_hist$mids);
#predictions = exp(predict(area_hist_model))
#lines(area_hist$mids,predictions, col='red', lwd=2)

mtext('D',adj=-.2,side=3,line=-0.25,cex=1.5)
hist(ind_exp_filt$ad_sig, main="", ylab = "FA Count", xlab = "Normalized Average Paxillin Intensity")

par(bty='n', mar=c(4,4.2,1.2,0.1))
hist(ind_exp_filt$ax[ind_exp_filt$ax < 8], main="", ylab = "FA Count",  xlab = "Axial Ratio")
hist(ind_exp_filt$cent_dist, main="", ylab = "FA Count",  
	 xlab = expression(paste("Distance from Edge (", symbol("m"), m, ')',sep='')))
hist(results_props$longevity, main="", ylab = "FA Count",  xlab = "Longevity (min)")
graphics.off()

svg(file.path(out_folder,'statics','longev_inset.svg'), width=3, height=3/2)
par(bty='n', mar=c(2,2,0.5,0))
hist(results_props$longevity[!is.na(results_props$longevity) & results_props$longevity > 20], 
    main="", ylab = "", xlab = "")
graphics.off()

print('Done with Static Properties')

########################################
#Kinetics Figure
#######################################
dir.create(dirname(file.path(out_folder,'kinetics','kinetics.svg')), 
    recursive=TRUE, showWarnings=FALSE);
svg(file.path(out_folder,'kinetics','kinetics.svg'),height=10.5);
layout(rbind(c(1,2),c(3,4),c(5,5)))

exp_one_only = load_results(exp_dirs[[1]],'intensity.Rdata')
exp_one_only = exp_one_only[[1]]

par(bty='n', mar=c(4,4.2,1.5,0))

plot.new()
mtext('A',adj=-.19,side=3,line=-0.5,cex=1.5)

ad_num = 1799
ad_num = 675
plot_ad_seq(exp_one_only, ad_num, type='overall', 
	phase_lengths=c(exp_one_only$assembly$length[ad_num],exp_one_only$disassembly$length[ad_num]))
mtext('B',adj=-.19,side=3,line=-0.5,cex=1.5)

par(bty='n', mar=c(4,4.2,4,0))
plot_ad_seq(exp_one_only, ad_num, main = 'Assembly');
limits = par("usr");
text(3,(limits[4]-limits[3])*0.8+limits[3],pos=3,expression(paste(R^2,' = 0.920')))
text(3,(limits[4]-limits[3])*0.8+limits[3],pos=3, offset=c(-0.65,0),
	substitute(paste('Slope = ', x), list(x=sprintf('%.03f',exp_one_only$assembly$slope[ad_num]))))
mtext('C',adj=-.19,side=3,line=-0.5,cex=1.5)

plot_ad_seq(exp_one_only,ad_num,type='disassembly', main = 'Disassembly')
limits = par("usr");
text(3,(limits[4]-limits[3])*0.8+limits[3],pos=3, expression(paste(R^2,' = 0.961')))
text(3,(limits[4]-limits[3])*0.8+limits[3], offset=c(-0.65,0),
	substitute(paste('Slope = ', x), list(x=sprintf('%.03f',exp_one_only$disassembly$slope[ad_num]))))
mtext('D',adj=-.19,side=3,line=-0.5,cex=1.5)
par(bty='n', mar=c(2.1,4.2,1.1,0))

boxplot_with_points(list(results_filt$a$slope,results_filt$dis$slope), 
		    names=c('Assembly', 'Disassembly'), boxwex=0.6, 
		    ylab=expression(paste('Rate (',min^-1,')',sep='')), point_cex=0.6)
#95% confidence intervals on the mean from Webb 2004
#segments(1.4,0.04,1.4,0.2,lwd=2)
#segments(1.35,0.12,1.45,0.12,lwd=2)
#segments(2.4,0.08,2.4,.088+0.004*2,lwd=2)
#segments(2.35,0.088,2.45,.088,lwd=2)
mtext('E',adj=-0.085,side=3,line=-0.5,cex=1.5)
graphics.off()

####################
#Supplemental
####################
dir.create(dirname(file.path(out_folder,'supplemental','R_squared.svg')), 
    recursive=TRUE, showWarnings=FALSE);
svg(file.path(out_folder,'supplemental','R_squared.svg'), width=14*.55, height=7*.55)
layout(rbind(c(1,2)))

par(bty='n', mar=c(4,4.2,2,0))

hist(results_nofilt$a$R, main='Assembly', freq=TRUE,
	 xlab=paste('Adjusted R Squared Values (n=',length(results_nofilt$a$R_sq),')', sep=''), 
	 ylab='# of Focal Adhesions')

plot_dims = par("usr");
sorted_r_vals = sort(results_nofilt$assem$R)
segments(sorted_r_vals[floor(length(sorted_r_vals)/2)],0, 
         sorted_r_vals[floor(length(sorted_r_vals)/2)],plot_dims[4], col='red', lwd=2)
mtext('A',adj=-.2,side=3,line=-0.5,cex=1.5)

hist(results_nofilt$dis$R,main='Disassembly', freq=TRUE,
	 xlab=paste('Adjusted R Squared Values (n=',length(results_nofilt$d$R_sq),')', sep=''), 
	 ylab='# of Focal Adhesions')

plot_dims = par("usr");
sorted_r_vals = sort(results_nofilt$dis$R)
segments(sorted_r_vals[floor(length(sorted_r_vals)/2)],0, 
         sorted_r_vals[floor(length(sorted_r_vals)/2)],plot_dims[4], col='red', lwd=2)
mtext('B',adj=-.2,side=3,line=-0.5,cex=1.5)
graphics.off()


print('Done with Kinetics')

########################################
#Spacial Figure
########################################

dir.create(dirname(file.path(out_folder,'spatial','spatial.svg')), 
    recursive=TRUE, showWarnings=FALSE);
svg(file.path(out_folder,'spatial','spatial.svg'))
par(bty='n',mar=c(4.2,4.1,0.1,1))
layout(rbind(c(1,2),c(3,4)))

breaks_end = ceil(max(c(results_filt$a$edge_dist,results_filt$dis$edge_dist), na.rm=T));
if (breaks_end %% 2 != 0) {
	breaks_end = breaks_end + 1;
}
these_breaks = seq(0,breaks_end,by=2);

hist(results_filt$a$edge_dist, 
     xlab=expression(paste('Distance from Edge at Birth (', mu, 'm) n=309', sep='')), 
     main = '', ylab = '# of Focal Adhesions', breaks=these_breaks)
mtext('A',adj=-.25,side=3,line=-1.5,cex=1.5)

par(bty='n',mar=c(4.2,4.1,0.1,0))
hist(results_filt$dis$edge_dist,
     xlab=expression(paste('Distance from Edge at Death (', mu, 'm) n=470', sep='')), 
     main = '', ylab = '# of Focal Adhesions', breaks=these_breaks)
mtext('B',adj=-.25,side=3,line=-1.5,cex=1.5)

par(bty='n',mar=c(4.2,4.1,0.1,1))
max_rate = max(c(results_filt$a$slope, results_filt$d$slope));
plot(results_filt$a$edge_dist, pch=19, cex=0.5,
	 results_filt$a$slope,
	 xlim = c(0,breaks_end),
         ylim = c(0,max_rate),
	 ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')),
	 xlab=expression(paste('Distance from Edge at Birth (', mu, 'm)', sep='')))
mtext('C',adj=-.25,side=3,line=-1.5,cex=1.5)

par(bty='n',mar=c(4.2,4.1,0.1,0))
plot(results_filt$d$edge_dist, pch=19, cex = 0.5,
	 results_filt$d$slope, 
	 xlim = c(0,breaks_end),
         ylim = c(0,max_rate),
	 ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')),
	 xlab=expression(paste('Distance from Edge at Death (', mu, 'm)', sep='')))
mtext('D',adj=-.25,side=3,line=-1.5,cex=1.5)
graphics.off()

plot(results_filt$a$edge_dist, pch=19, cex=0.5,
	 results_filt$a$slope,
	 xlim = c(0,breaks_end),
         ylim = c(0,max_rate),
	 ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')),
	 xlab=expression(paste('Distance from Edge at Birth (', mu, 'm)', sep='')))
axis_ticks = axTicks(2);
graphics.off();

#Alternate combined plots
dir.create(dirname(file.path(out_folder,'spatial','spatial_alt.svg')), 
    recursive=TRUE, showWarnings=FALSE);
svg(file.path(out_folder,'spatial','spatial_alt.svg'),height=4, width=10)
par(bty='n',mar=c(4.2,4.1,0.1,6))
layout(rbind(c(1,2)))

breaks_end = ceil(max(c(results_filt$a$edge_dist,results_filt$dis$edge_dist), na.rm=T));
if (breaks_end %% 2 != 0) {
	breaks_end = breaks_end + 1;
}
these_breaks = seq(0,breaks_end,by=2);
max_rate = max(c(results_filt$a$slope, results_filt$d$slope));

#assembly rates
hist(results_filt$a$edge_dist, 
     xlab=expression(paste('Distance from Edge at Birth (', mu, 'm) n=309', sep='')), 
     main = '', ylab = '# of Focal Adhesions', breaks=these_breaks)
plot_range = par('usr')
# segments(1.5,0,1.5,plot_range[4], col='purple', lwd=2)
points(results_filt$a$edge_dist, pch=19, cex=0.35, col='darkgreen',
	 results_filt$a$slope*((plot_range[4]-10)/max_rate),
	 ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')))
# axis_ticks = axTicks(2);
axis(4, at = axis_ticks*((plot_range[4]-10)/max_rate), labels=axis_ticks, col='darkgreen')
mtext(expression(paste('Assembly Rate (',min^-1,')',sep='')),side=4,line=3, col='darkgreen');
mtext('A',adj=-.22,side=3,line=-1.5,cex=1.75)

#disassembly rates
par(bty='n',mar=c(4.2,4.1,0.1,4))
hist(results_filt$dis$edge_dist,
     xlab=expression(paste('Distance from Edge at Death (', mu, 'm) n=470', sep='')), 
     main = '', ylab = '# of Focal Adhesions', breaks=these_breaks)
plot_range = par('usr')
# segments(1.5,0,1.5,plot_range[4], col='purple', lwd=2)
points(results_filt$d$edge_dist, pch=19, cex=0.35, col='red',
	 results_filt$d$slope*((plot_range[4]-10)/max_rate),
	 ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')))
axis(4, at = axis_ticks*((plot_range[4]-10)/max_rate), labels=axis_ticks, col='red')
mtext(expression(paste('Disassembly Rate (',min^-1,')',sep='')),side=4,line=3, col='red');
mtext('B',adj=-0.22,side=3,line=-1.5,cex=1.75)

graphics.off()

####################
#Supplemental
####################
dir.create(dirname(file.path(out_folder,'supplemental','spatial_nofilt.svg')), 
    recursive=TRUE, showWarnings=FALSE);
svg(file.path(out_folder,'supplemental','spatial_nofilt.svg'))
par(bty='n',mar=c(4.2,4.1,2,0))
layout(rbind(c(1,2),c(3,4)))

breaks_end = ceil(max(c(results_nofilt$a$edge_dist,results_nofilt$dis$edge_dist), na.rm=T));
if (breaks_end %% 2 != 0) {
	breaks_end = breaks_end + 1;
}
these_breaks = seq(0,breaks_end,by=2);

hist(results_nofilt$a$edge_dist, 
    xlab=expression(paste('Distance from Edge at Birth (', mu, 'm) n=1169', sep='')),
    main = '', ylab = '# of Focal Adhesions', breaks=these_breaks)
mtext('A',adj=-.2,side=3,line=-0.5,cex=1.5)

hist(results_nofilt$dis$edge_dist,
     xlab=expression(paste('Distance from Edge at Death (', mu, 'm) n=1460', sep='')), 
     main = '', ylab = '# of Focal Adhesions',breaks=these_breaks)
mtext('B',adj=-.2,side=3,line=-0.5,cex=1.5)

par(bty='n',mar=c(4.2,4.1,2,0.4))

smoothScatter(results_nofilt$a$edge_dist, 
	 results_nofilt$a$slope,
	 xlim = c(0,breaks_end), 
	 ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')),
	 xlab=expression(paste('Distance from Edge at Birth (', mu, 'm)', sep='')))
mtext('C',adj=-.2,side=3,line=-0.5,cex=1.5)

par(bty='n',mar=c(4.2,4.1,2,0.95))

smoothScatter(results_nofilt$d$edge_dist, 
	 results_nofilt$d$slope,
	 xlim = c(0,breaks_end), 
	 ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')),
	 xlab=expression(paste('Distance from Edge at Death (', mu, 'm)', sep='')))
mtext('D',adj=-.2,side=3,line=-0.5,cex=1.5)
graphics.off()


svg(file.path(out_folder,'supplemental','birth_vs_death_pos.svg'))
#pdf(file.path(out_folder,'supplemental','birth_vs_death_pos.pdf'))
par(bty='n',mar=c(4.2,4.1,2,0.2))
plot(results_nofilt$j$birth_dist, results_nofilt$j$death_dist, 
     xlab=expression(paste('Distance from Edge at Birth (', mu, 'm)', sep='')),
     ylab=expression(paste('Distance from Edge at Death (', mu, 'm)', sep='')), 
     pch=20, cex=0.75)
#points(results_S_nofilt$j$b, results_S_nofilt$j$d, col='red', pch=20, cex=0.5)

birth_vs_death_model <- lm(death_dist ~ birth_dist, data=results_nofilt$joint)
abline(birth_vs_death_model, col='green', lwd = 3)
model_summary <- summary(birth_vs_death_model)
x_data <- data.frame(birth_dist = seq(min(results_nofilt$j$b, na.rm=T),max(results_nofilt$j$b, na.rm=T),by=0.01))
line_conf = predict(birth_vs_death_model, x_data, interval="confidence", level=0.95)
lines(x_data$birth_dist, line_conf[,2], col='red', lty=2, lwd = 3)
lines(x_data$birth_dist, line_conf[,3], col='red', lty=2, lwd = 3)

segments(0,0,max(results_nofilt$j$b, na.rm=T), max(results_nofilt$j$b, na.rm=T), col='blue', lty=4, lwd = 3)

graphics.off()

svg(file.path(out_folder,'supplemental','birth_vs_death_pos_S178A.svg'))
#pdf(file.path(out_folder,'supplemental','birth_vs_death_pos.pdf'))
par(bty='n',mar=c(4.2,4.1,2,0.2))
plot(results_S_nofilt$j$birth_dist, results_S_nofilt$j$death_dist, 
     xlab=expression(paste('Distance from Edge at Birth (', mu, 'm)', sep='')),
     ylab=expression(paste('Distance from Edge at Death (', mu, 'm)', sep='')), 
     pch=20, cex=0.75)
#points(results_S_nofilt$j$b, results_S_nofilt$j$d, col='red', pch=20, cex=0.5)

birth_vs_death_model <- lm(death_dist ~ birth_dist, data=results_S_nofilt$joint)
abline(birth_vs_death_model, col='green', lwd = 3)
model_summary <- summary(birth_vs_death_model)
x_data <- data.frame(birth_dist = seq(min(results_S_nofilt$j$b, na.rm=T),max(results_S_nofilt$j$b, na.rm=T),by=0.01))
line_conf = predict(birth_vs_death_model, x_data, interval="confidence", level=0.95)
lines(x_data$birth_dist, line_conf[,2], col='red', lty=2, lwd = 3)
lines(x_data$birth_dist, line_conf[,3], col='red', lty=2, lwd = 3)

segments(0,0,max(results_S_nofilt$j$b, na.rm=T), max(results_S_nofilt$j$b, na.rm=T), col='blue', lty=4, lwd = 3)

graphics.off()
print('Done with Spacial')

############################################################
#Comparing S178A to Wild-type
############################################################
dir.create(dirname(file.path(out_folder,'lifetimes','stable_mean.pdf')), 
    recursive=TRUE, showWarnings=FALSE);
pdf(file.path(out_folder,'lifetimes','stable_mean.pdf'), width=8.5, pointsize=8)
boxplot_with_points(
        list(corr_results_onlysignif$joint$stable_mean, 
            corr_results_S_onlysignif$joint$stable_mean, 
            #CB_results_onlysignif$joint$stable_mean, 
            #CB_results_S_onlysignif$joint$stable_mean,
            results_onlysignif$joint$stable_mean, 
            results_S_onlysignif$joint$stable_mean
            ),
        #names=c('Local WT','Local S189A','CB WT', 'CB S178A', 'WT','S178A'), ylab='Mean Stable Intensity')
        names=c('Local WT','Local S189A', 'WT','S178A'), ylab='Mean Stable Intensity')
graphics.off()

print('Done with Stable Lifetime Averages')

########################################
#Lifetime Phases
########################################

stage_data <- gather_stage_lengths(corr_results_onlysignif, corr_results_S_onlysignif)

dir.create(dirname(file.path(out_folder,'lifetimes','adhesion_phase_lifetimes.svg')), 
    recursive=TRUE, showWarnings=FALSE);
svg(file.path(out_folder,'lifetimes','adhesion_phase_lifetimes.svg'))
par(bty='n',mar=c(2,4,0,0))

err_bars = plot_stage_length_data(stage_data, top_gap=3.5, names=c('Wild-type','S178A'))

#Signifcance Bars
bar_length = 1;
sep_from_data = 0.5;

#upside down identification bar
upper_left = c(err_bars[1,1], min(err_bars[1,3],err_bars[4,3]) - sep_from_data);
lower_right = c(err_bars[4,1], min(err_bars[1,3],err_bars[4,3]) - (sep_from_data + bar_length));
lines(c(upper_left[1],upper_left[1],lower_right[1], lower_right[1]),
	  c(upper_left[2],lower_right[2],lower_right[2],upper_left[2]))  
text(mean(c(upper_left[1],lower_right[1]))-0.005,lower_right[2]-sep_from_data,"**",cex=1.5)

#right side up identification bar
upper_left = c(err_bars[3,1], max(err_bars[3,4],err_bars[6,4]) + sep_from_data + bar_length);
lower_right = c(err_bars[6,1], max(err_bars[3,4],err_bars[6,4]) + sep_from_data);
lines(c(upper_left[1],upper_left[1],lower_right[1], lower_right[1]),
	  c(lower_right[2],upper_left[2],upper_left[2],lower_right[2]))  
text(mean(c(upper_left[1],lower_right[1]))-0.005,upper_left[2]+sep_from_data,"*",cex=1.5)

graphics.off()

svg(file.path(out_folder,'lifetimes','adhesion_phase_lifetimes_alt.svg'))
par(bty='n',mar=c(2.1,4,0,0))

sideways_err_bars = plot_stage_length_data(stage_data, type='side_by_side', top_gap=1, 
    names=c('Wild-type','S178A'), sideways_high_xlim=11);

#Signifcance Bars
bar_length = 0.4;
sep_from_data = 0.2;

#assembly significance bar
upper_left = c(sideways_err_bars[1,1], max(sideways_err_bars[1,4],sideways_err_bars[4,4]) + sep_from_data + bar_length);
lower_right = c(sideways_err_bars[4,1], max(sideways_err_bars[1,4],sideways_err_bars[4,4]) + sep_from_data);
lines(c(upper_left[1],upper_left[1],lower_right[1], lower_right[1]),
	  c(lower_right[2],upper_left[2],upper_left[2],lower_right[2]))  
text(mean(c(upper_left[1],lower_right[1]))-0.005,upper_left[2]+0.2,"**",cex=1.5)

#disassembly significance bar
upper_left = c(sideways_err_bars[3,1], max(sideways_err_bars[3,4],sideways_err_bars[6,4]) + sep_from_data + bar_length);
lower_right = c(sideways_err_bars[6,1], max(sideways_err_bars[3,4],sideways_err_bars[6,4]) + sep_from_data);
lines(c(upper_left[1],upper_left[1],lower_right[1], lower_right[1]),
	  c(lower_right[2],upper_left[2],upper_left[2],lower_right[2]))  
text(mean(c(upper_left[1],lower_right[1]))-0.005,upper_left[2]+sep_from_data,"*",cex=1.5)

graphics.off();

print('Done with Lifetime Phase Lengths')

########################################
#Dynamics
########################################
#svg(file.path(out_folder,'S178A','S178A_vs_wild-type.svg'))
dir.create(dirname(file.path(out_folder,'S178A','S178A_vs_wild-type.pdf')), 
    recursive=TRUE, showWarnings=FALSE);
pdf(file.path(out_folder,'S178A','S178A_vs_wild-type.pdf'))
layout(rbind(c(1,2),c(3,4)))
par(bty='n', mar=c(2,4,0,0))

max_rate = max(c(results_filt$as$slope,results_S_filt$as$slope,
                 results_filt$dis$slope,results_S_filt$dis$slope));

#Panel Assembly Rates
boxplot_with_points(list(results_filt$as$slope,results_S_filt$as$slope), 
        names=c('Wild-type','S178A'), 
        colors=c('orange','blue'),
        ylim = c(0,max_rate + 0.012),
        ylab=expression(paste('Assembly Rate (',min^-1,')',sep=''))
)

bar_length = .005;
sep_from_data = 0.005;

upper_left = c(1, max_rate + sep_from_data + bar_length);
lower_right = c(2, max_rate + sep_from_data);
lines(c(upper_left[1],upper_left[1],lower_right[1], lower_right[1]),
	  c(lower_right[2],upper_left[2],upper_left[2],lower_right[2]))  
text(mean(c(upper_left[1],lower_right[1]))-0.005,upper_left[2]+sep_from_data,"***",cex=1.5)
mtext('A',adj=-.25,side=3,line=-1.5,cex=1.5)	    

#Panel Disassembly Rates
boxplot_with_points(list(results_filt$dis$slope,results_S_filt$dis$slope), 
    names=c('Wild-type','S178A'), 
    colors=c('orange','blue'),
    ylim = c(0,max_rate + 0.012),
    ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')) 
)

upper_left = c(1, max_rate + sep_from_data + bar_length);
lower_right = c(2, max_rate + sep_from_data);
lines(c(upper_left[1],upper_left[1],lower_right[1], lower_right[1]),
	  c(lower_right[2],upper_left[2],upper_left[2],lower_right[2]))  
text(mean(c(upper_left[1],lower_right[1]))-0.005,upper_left[2]+sep_from_data,"*",cex=1.5)
mtext('B',adj=-.25,side=3,line=-1.5,cex=1.5)

#Panel Birth Distances
par(bty='n', mar=c(2,4,1.5,0))
max_dist = max(c(results_filt$as$edge_dist,results_S_filt$as$edge_dist,
                 results_filt$dis$edge_dist,results_S_filt$dis$edge_dist), na.rm=T)
boxplot_with_points(list(results_filt$as$edge_dist,results_S_filt$as$edge_dist), 
    names=c('Wild-type','S178A'), 
    ylim=c(0,max_dist+2), 
    colors=c('orange','blue'),
    ylab=expression(paste('Distance from Edge at Birth (',mu,'m)',sep=''))
)

bar_length = 1;
sep_from_data = 1;

# upper_left = c(1, max_dist + sep_from_data + bar_length);
# lower_right = c(2, max_dist + sep_from_data);
# lines(c(upper_left[1],upper_left[1],lower_right[1], lower_right[1]),
# 	  c(lower_right[2],upper_left[2],upper_left[2],lower_right[2]))  
# text(mean(c(upper_left[1],lower_right[1]))-0.005,upper_left[2]+sep_from_data,"**",cex=1.5)
mtext('C',adj=-.25,side=3,line=-1.5,cex=1.5)

#Panel Death Distances
boxplot_with_points(list(results_filt$dis$edge_dist,results_S_filt$dis$edge_dist), 
    names=c('Wild-type','S178A'), 
    ylim=c(0,max_dist+2),
    colors=c('orange','blue'),
    ylab=expression(paste('Distance from Edge at Death (',mu,'m)',sep=''))
)
mtext('D',adj=-.25,side=3,line=-1.5,cex=1.5)	    

graphics.off()

#comparing S178A area results
dir.create(dirname(file.path(out_folder,'S178A','S178A_vs_wild-type.pdf')), 
    recursive=TRUE, showWarnings=FALSE);
pdf(file.path(out_folder,'S178A','S178A_vs_wild-type_area.pdf'))
layout(rbind(c(1,2),c(3,4)))
par(bty='n', mar=c(2,4,0,0))

max_rate = max(c(area_onlysignif$as$slope,area_S_onlysignif$as$slope,
                 area_onlysignif$dis$slope,area_S_onlysignif$dis$slope));

#Panel Assembly Rates
boxplot_with_points(list(area_onlysignif$as$slope,area_S_onlysignif$as$slope), 
        names=c('Wild-type','S178A'), 
        colors=c('orange','blue'),
        ylim = c(0,max_rate + 0.012),
        ylab=expression(paste('Assembly Rate (',min^-1,')',sep=''))
)

bar_length = .005;
sep_from_data = 0.005;

upper_left = c(1, max_rate + sep_from_data + bar_length);
lower_right = c(2, max_rate + sep_from_data);
lines(c(upper_left[1],upper_left[1],lower_right[1], lower_right[1]),
	  c(lower_right[2],upper_left[2],upper_left[2],lower_right[2]))  
text(mean(c(upper_left[1],lower_right[1]))-0.005,upper_left[2]+sep_from_data,"***",cex=1.5)
mtext('A',adj=-.25,side=3,line=-1.5,cex=1.5)	    

#Panel Disassembly Rates
boxplot_with_points(list(area_onlysignif$dis$slope,area_S_onlysignif$dis$slope), 
    names=c('Wild-type','S178A'), 
    colors=c('orange','blue'),
    ylim = c(0,max_rate + 0.012),
    ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')) 
)

upper_left = c(1, max_rate + sep_from_data + bar_length);
lower_right = c(2, max_rate + sep_from_data);
lines(c(upper_left[1],upper_left[1],lower_right[1], lower_right[1]),
	  c(lower_right[2],upper_left[2],upper_left[2],lower_right[2]))  
text(mean(c(upper_left[1],lower_right[1]))-0.005,upper_left[2]+sep_from_data,"*",cex=1.5)
mtext('B',adj=-.25,side=3,line=-1.5,cex=1.5)

#Panel Birth Distances
par(bty='n', mar=c(2,4,1.5,0))
max_dist = max(c(area_onlysignif$as$edge_dist,area_S_onlysignif$as$edge_dist,
                 area_onlysignif$dis$edge_dist,area_S_onlysignif$dis$edge_dist), na.rm=T)
boxplot_with_points(list(area_onlysignif$as$edge_dist,area_S_onlysignif$as$edge_dist), 
    names=c('Wild-type','S178A'), 
    ylim=c(0,max_dist+2), 
    colors=c('orange','blue'),
    ylab=expression(paste('Distance from Edge at Birth (',mu,'m)',sep=''))
)

bar_length = 1;
sep_from_data = 1;

upper_left = c(1, max_dist + sep_from_data + bar_length);
lower_right = c(2, max_dist + sep_from_data);
lines(c(upper_left[1],upper_left[1],lower_right[1], lower_right[1]),
	  c(lower_right[2],upper_left[2],upper_left[2],lower_right[2]))  
text(mean(c(upper_left[1],lower_right[1]))-0.005,upper_left[2]+sep_from_data,"**",cex=1.5)
mtext('C',adj=-.25,side=3,line=-1.5,cex=1.5)

#Panel Death Distances
boxplot_with_points(list(area_onlysignif$dis$edge_dist,area_S_onlysignif$dis$edge_dist), 
    names=c('Wild-type','S178A'), 
    ylim=c(0,max_dist+2),
    colors=c('orange','blue'),
    ylab=expression(paste('Distance from Edge at Death (',mu,'m)',sep=''))
)
mtext('D',adj=-.25,side=3,line=-1.5,cex=1.5)	    

graphics.off()

###############
#Conf Levels
###############
require(boot)

determine_median_p_value(results_S_filt$as$slope, results_filt$as$slope)

determine_median_p_value(results_S_filt$dis$slope, results_filt$dis$slope)

determine_median_p_value(results_filt$as$edge_dist, results_S_filt$as$edge_dist)

determine_median_p_value(results_filt$dis$edge_dist, results_S_filt$dis$edge_dist)

###############
#Supplemental
###############
dir.create(dirname(file.path(out_folder,'supplemental','unfilt_S178A_vs_wild-type.pdf')), 
    recursive=TRUE, showWarnings=FALSE);
pdf(file.path(out_folder,'supplemental','unfilt_S178A_vs_wild-type.pdf'))
layout(rbind(c(1,2),c(3,4)))
par(bty='n', mar=c(2,4,0,0))

max_rate = max(c(results_nofilt$as$slope,results_S_nofilt$as$slope,
                 results_nofilt$dis$slope,results_S_nofilt$dis$slope));

#Panel Assembly Rates
boxplot_with_points(list(results_nofilt$as$slope,results_S_nofilt$as$slope), 
        names=c('Wild-type','S178A'), 
        colors=c('orange','blue'),
        ylim = c(0,max_rate),
        ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')),
        inc.points=FALSE
)

#bar_length = .005;
#sep_from_data = 0.005;
#
#upper_left = c(1, max_rate + sep_from_data + bar_length);
#lower_right = c(2, max_rate + sep_from_data);
#lines(c(upper_left[1],upper_left[1],lower_right[1], lower_right[1]),
#	  c(lower_right[2],upper_left[2],upper_left[2],lower_right[2]))  
#text(mean(c(upper_left[1],lower_right[1]))-0.005,upper_left[2]+sep_from_data,"***",cex=1.5)
mtext('A',adj=-.25,side=3,line=-1.5,cex=1.5)	    

#Panel Disassembly Rates
boxplot_with_points(list(results_nofilt$dis$slope,results_S_nofilt$dis$slope), 
    names=c('Wild-type','S178A'), 
    colors=c('orange','blue'),
    ylim = c(0,max_rate),
    ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')),
    inc.points = FALSE
)

#upper_left = c(1, max_rate + sep_from_data + bar_length);
#lower_right = c(2, max_rate + sep_from_data);
#lines(c(upper_left[1],upper_left[1],lower_right[1], lower_right[1]),
#	  c(lower_right[2],upper_left[2],upper_left[2],lower_right[2]))  
#text(mean(c(upper_left[1],lower_right[1]))-0.005,upper_left[2]+sep_from_data,"*",cex=1.5)
mtext('B',adj=-.25,side=3,line=-1.5,cex=1.5)

#Panel Birth Distances
par(bty='n', mar=c(2,4,1.5,0))
max_dist = max(c(results_nofilt$as$edge_dist,results_S_nofilt$as$edge_dist,
                 results_nofilt$dis$edge_dist,results_S_nofilt$dis$edge_dist), na.rm=T)
boxplot_with_points(list(results_nofilt$as$edge_dist,results_S_nofilt$as$edge_dist), 
    names=c('Wild-type','S178A'), 
    ylim=c(0,max_dist), 
    colors=c('orange','blue'),
    ylab=expression(paste('Distance from Edge at Birth (',mu,'m)',sep='')),
    inc.points = FALSE
)

#bar_length = 1;
#sep_from_data = 1;
#
#upper_left = c(1, max_dist + sep_from_data + bar_length);
#lower_right = c(2, max_dist + sep_from_data);
#lines(c(upper_left[1],upper_left[1],lower_right[1], lower_right[1]),
#	  c(lower_right[2],upper_left[2],upper_left[2],lower_right[2]))  
#text(mean(c(upper_left[1],lower_right[1]))-0.005,upper_left[2]+sep_from_data,"**",cex=1.5)
mtext('C',adj=-.25,side=3,line=-1.5,cex=1.5)

#Panel Death Distances
boxplot_with_points(list(results_nofilt$dis$edge_dist,results_S_nofilt$dis$edge_dist), 
    names=c('Wild-type','S178A'), 
    ylim=c(0,max_dist),
    colors=c('orange','blue'),
    ylab=expression(paste('Distance from Edge at Death (',mu,'m)',sep='')), 
    inc.points = FALSE
)
mtext('D',adj=-.25,side=3,line=-1.5,cex=1.5)	    

graphics.off()
print('Done with S178A Comparisons')

############################################################
#Distance versus Background corrected correlation
############################################################

all = list()
for (i in 1:length(background_correlation)) {
	all$distances = c(all$distances, background_correlation[[i]]$distances)
	all$correlations = c(all$correlations, background_correlation[[i]]$correlations)
}	
all = as.data.frame(all)

all_S = list()
for (i in 1:length(background_correlation_S)) {
	all_S$distances = c(all_S$distances, background_correlation_S[[i]]$distances)
	all_S$correlations = c(all_S$correlations, background_correlation_S[[i]]$correlations)
}	
all_S = as.data.frame(all_S)

binned_all = bin_corr_data(all, bin_max=200)
binned_all_S = bin_corr_data(all_S, bin_max=200)

dir.create(dirname(file.path(out_folder,'S178A','dist_vs_corr.pdf')), 
    recursive=TRUE, showWarnings=FALSE);
pdf(file.path(out_folder,'S178A','dist_vs_corr.pdf'))
y_max = max(c(binned_all$upper, binned_all_S$upper), na.rm=TRUE)

plot(binned_all$mids, binned_all$means, ylim=c(0, y_max), 
        xlab=expression(paste('Mean Distance (', mu, 'm)', sep='')), ylab="correlation")
errbar(binned_all$mids, binned_all$means, binned_all$upper, binned_all$lower, add=TRUE, xlab="", ylab="")

#plot(binned_all_S$mids, binned_all_S$means, ylim=c(0, y_max), col='red', add=TRUE)
errbar(binned_all_S$mids, binned_all_S$means, binned_all_S$upper, binned_all_S$lower, add=TRUE,col='red', xlab="", ylab="")
graphics.off()

print('Done with Distance versus Pax Concentration Correlation')
