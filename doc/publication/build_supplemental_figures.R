rm(list = ls())
source('FA_analysis_lib.R')
library(lattice)
library(geneplotter)
library(Hmisc)

################################################################################
#Result loading
################################################################################
exp_dirs <- Sys.glob('../../results/focal_adhesions/*/adhesion_props/models/')
#exp_dirs <- Sys.glob('../../results/lin_region_variation/FA/6/*/adhesion_props/models/')
#exp_dirs <- Sys.glob('../../results/lin_region_variation/FA/8/*/adhesion_props/models/')
#exp_dirs <- Sys.glob('../../results/lin_region_variation/FA/no_log_trans/*/adhesion_props/models/')
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
results = load_results(exp_dirs,file.path('intensity_model.Rdata'))
area = load_results(exp_dirs,file.path('area_model.Rdata'))
ind_results <- load_data_files(exp_dirs,file.path('..','individual_adhesions.csv'), headers=T, debug=FALSE, inc_exp_names=FALSE);

cell_area <- load_data_files(exp_dirs,file.path('..','Cell_size.csv'), headers=F, debug=FALSE, inc_exp_names=FALSE);
ad_area <- load_data_files(exp_dirs,file.path('..','Ad_size.csv'), headers=F, debug=FALSE, inc_exp_names=FALSE);

exp_dirs_S <- Sys.glob('../../results/S178A/*/adhesion_props/models/')
#exp_dirs_S <- Sys.glob('../../results/lin_region_variation/S178A/6/*/adhesion_props/models/')
#exp_dirs_S <- Sys.glob('../../results/lin_region_variation/S178A/8/*/adhesion_props/models/')
exp_dirs_S <- exp_dirs_S[file_test('-d',exp_dirs_S)]
results_S = load_results(exp_dirs_S,file.path('intensity_model.Rdata'))
area_S = load_results(exp_dirs_S,file.path('area_model.Rdata'))
ind_results_S <- load_data_files(exp_dirs_S,file.path('..','individual_adhesions.csv'), headers=T, debug=FALSE, inc_exp_names=FALSE);

cell_area_S <- load_data_files(exp_dirs_S,file.path('..','Cell_size.csv'), headers=F, debug=FALSE, inc_exp_names=FALSE);
ad_area_S <- load_data_files(exp_dirs_S,file.path('..','Ad_size.csv'), headers=F, debug=FALSE, inc_exp_names=FALSE);

print('Done Loading Data')

########################################
#Result filtering
########################################
results_props = gather_general_props(results)
ind_exp_filt = gather_single_image_props(ind_results)
results_nofilt = filter_results(results, min_R_sq = -Inf, max_p_val = Inf)
results_S_nofilt = filter_results(results_S, min_R_sq = -Inf, max_p_val = Inf)
area_nofilt = filter_mixed_area(area, results, min_R_sq = -Inf, max_p_val = Inf)

results_S_props = gather_general_props(results_S)
ind_exp_filt_S = gather_single_image_props(ind_results_S)
results_filt = filter_results(results)
results_S_filt = filter_results(results_S)
area_filt = filter_mixed_area(area, results)
print('Done Filtering Data')

################################################################################
#Plotting
################################################################################
out_folder = '../../doc/publication/figures/'

########################################
#Comparing Area/Kinetics
########################################
svg(file.path(out_folder,'supplemental','area_vs_intensity_rsq.svg'),height=12)
layout(rbind(c(1,2),c(3),c(4)), height=c(1,2,2))

par(bty='n', mar=c(5,4.2,2,0.1))
hist(area_nofilt$a$R, main="Assembling Area R squared", xlab='')
mtext('A',adj=0,cex=1.5)
hist(area_nofilt$d$R, main="Dissembling Area R squared", xlab='')
mtext('B',adj=0,cex=1.5)

plot(results_nofilt$a$R_sq,area_nofilt$a$R, xlab="Assembling Intensity R squared", ylab = "Assembling Area R squared", main='Pearson Correlation (95% Conf): 0.53-0.59')
mtext('C',adj=0,cex=1.5)

plot(results_nofilt$d$R_sq,area_nofilt$d$R, xlab="Disassembling Intensity R squared", ylab = "Disassembling Area R squared", main='Pearson Correlation (95% Conf): 0.58-0.64')
mtext('D',adj=0,cex=1.5)
graphics.off()

svg(file.path(out_folder,'area_vs_intensity_slope.svg'))
layout(rbind(c(1,2),c(3,4)))
par(bty='n', mar=c(5,4.2,2,1))

smoothScatter(results_nofilt$a$sl,area_nofilt$a$sl, nbin=64, xlab='Intensity Slope', ylab='Area Slope', main='Assembly')
abline(lm(area_nofilt$a$sl~results_nofilt$a$sl), col='red')
mtext('A',adj=0,cex=1.5)

plot(results_filt$a$sl,area_filt$a$sl, xlab='Filtered Intensity Slope', ylab='Filtered Area Slope', main='Assembly')
abline(lm(area_filt$a$sl~results_filt$a$sl), col='red')
mtext('B',adj=0,cex=1.5)

smoothScatter(results_nofilt$d$sl,area_nofilt$d$sl, nbin=64, xlab='Intensity Slope', ylab='Area Slope', main='Disassembly')
abline(lm(area_nofilt$d$sl~results_nofilt$d$sl), col='red')
mtext('C',adj=0,cex=1.5)

plot(results_filt$d$sl,area_filt$d$sl, xlab='Filtered Intensity Slope', ylab='Filtered Area Slope', main ='Disassembly')
abline(lm(area_filt$d$sl~results_filt$d$sl),col='red')
mtext('D',adj=0,cex=1.5)
graphics.off()
print('Done with Area Properties')

########################################
#Comparing S178A to Wild-type
########################################
svg(file.path(out_folder,'supplemental','unfilt_S178A_vs_wild-type.svg'))

layout(rbind(c(1,2),c(3,4)))
par(bty='n', mar=c(5,4.2,2,0.1))

boxplot_with_points(list(results_nofilt$as$slope,results_S_nofilt$as$slope), inc.points = FALSE,
				    names=c('Wild-type','S178A'),					ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')))
mtext('A',adj=0,cex=1.5)				    

boxplot_with_points(list(results_nofilt$dis$slope,results_S_nofilt$dis$slope), inc.points = FALSE,
				    names=c('Wild-type','S178A'), 
				    ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')))
mtext('B',adj=0,cex=1.5)

boxplot_with_points(list(results_nofilt$as$edge_dist,results_S_nofilt$as$edge_dist), inc.points = FALSE,
					names=c('Wild-type','S178A'), 
					ylab=expression(paste('Distance from Edge at Birth (',mu,'m)',sep='')))
mtext('C',adj=0,cex=1.5)

boxplot_with_points(list(results_nofilt$dis$edge_dist,results_S_nofilt$dis$edge_dist), inc.points = FALSE,
					names=c('Wild-type','S178A'), 
					ylab=expression(paste('Distance from Edge at Death (',mu,'m)',sep='')))
mtext('D',adj=0,cex=1.5)					

graphics.off()
print('Done Mutant Comparison Figures')