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
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
results = load_results(exp_dirs,file.path('intensity_model.Rdata'))
area = load_results(exp_dirs,file.path('area_model.Rdata'))

ind_results <- load_data_files(exp_dirs,file.path('..','individual_adhesions.csv'), T, debug=FALSE, inc_exp_names=FALSE);

exp_dirs_S <- Sys.glob('../../results/S178A/*/adhesion_props/models/')
#exp_dirs_S <- Sys.glob('../../results/lin_region_variation/S178A/6/*/adhesion_props/models/')
#exp_dirs_S <- Sys.glob('../../results/lin_region_variation/S178A/8/*/adhesion_props/models/')
exp_dirs_S <- exp_dirs_S[file_test('-d',exp_dirs_S)]
results_S = load_results(exp_dirs_S,file.path('intensity_model.Rdata'))
area_S = load_results(exp_dirs_S,file.path('area_model.Rdata'))

ind_results_S <- load_data_files(exp_dirs_S,file.path('..','individual_adhesions.csv'), T, debug=FALSE, inc_exp_names=FALSE);
print('Done Loading Data')

########################################
#Result filtering
########################################
results_props = gather_general_props(results)
ind_exp_filt = gather_single_image_props(ind_results)
results_nofilt = filter_mixed_results(results, results, min_R_sq = -Inf, max_p_val = Inf)
results_S_nofilt = filter_mixed_results(results_S, results_S, min_R_sq = -Inf, max_p_val = Inf)
area_nofilt = filter_mixed_area(area, results, min_R_sq = -Inf, max_p_val = Inf)

results_S_props = gather_general_props(results_S)
ind_exp_filt_S = gather_single_image_props(ind_results_S)
results_filt = filter_mixed_results(results, results)
results_S_filt = filter_mixed_results(results_S, results_S)
area_filt = filter_mixed_area(area, results)
print('Done Filtering Data')

################################################################################
#Plotting
################################################################################
out_folder = '../../doc/publication/figures/'

########################################
#Comparing S178A to Wild-type
########################################
pdf(file.path(out_folder,'S178A_vs_wild-type.pdf'),height=9)

layout(rbind(c(1,2),c(3,4),c(5,6)))
par(bty='n', mar=c(5,4.2,2,0.1))

boxplot_with_points(list(ind_exp_filt$Area,ind_exp_filt_S$Area), inc.points = FALSE,
				    names=c('Wild-type','S178A'),
				    ylab = expression(paste('Adhesion Area (', symbol("m"), m^2, ')',sep='')))
mtext('A',adj=0,cex=1.5)

boxplot_with_points(list(ind_exp_filt$axial,ind_exp_filt_S$axial), inc.points = FALSE,
				    names=c('Wild-type','S178A'),
				    ylab = 'Axial Ratio')
mtext('B',adj=0,cex=1.5)

boxplot_with_points(list(results_filt$as$slope,results_S_filt$as$slope), 
				    names=c('Wild-type','S178A'),					ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')))
mtext('C',adj=0,cex=1.5)				    

boxplot_with_points(list(results_filt$dis$slope,results_S_filt$dis$slope), 
				    names=c('Wild-type','S178A'), 
				    ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')))
mtext('D',adj=0,cex=1.5)

boxplot_with_points(list(results_filt$as$edge_dist,results_S_filt$as$edge_dist), 
					names=c('Wild-type','S178A'), 
					ylab=expression(paste('Distance from Edge at Birth (',mu,'m)',sep='')))
mtext('E',adj=0,cex=1.5)					

boxplot_with_points(list(results_filt$dis$edge_dist,results_S_filt$dis$edge_dist), 
					names=c('Wild-type','S178A'), 
					ylab=expression(paste('Distance from Edge at Death (',mu,'m)',sep='')))
mtext('F',adj=0,cex=1.5)					

graphics.off()

pdf(file.path(out_folder,'unfilt_S178A_vs_wild-type.pdf'))

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
stopifnot(FALSE)

########################################
#All by all plots
########################################
#pdf(file.path(out_folder,'as_splom.pdf'),width=12,height=12)
#splom(cbind(results_filt$a[1:10]), sub='Assembly Scatter Plot (Filtered)')
#splom(cbind(results_nofilt$a[1:10]), sub='Assembly Scatter Plot')
#graphics.off()
#
#pdf(file.path(out_folder,'dis_splom.pdf'),width=12,height=12)
#splom(cbind(results_filt$dis[1:10]), sub='Disassembly Scatter Plot (Filtered)')
#splom(cbind(results_nofilt$dis[1:10]), sub='Disassembly Scatter Plot')
#graphics.off()

########################################
#General Properties
########################################

pdf(file.path(out_folder,'general_props.pdf'))
layout(cbind(c(1,1,2,2,3,3), c(0,4,4,5,5,0)))
par(bty='n', mar=c(5,4.2,2,0.1))

hist(ind_exp_filt$Area, main="", ylab = "Adhesion Count", xlab = expression(paste('Adhesion Area (', symbol("m"), m^2, ')',sep='')))
mtext('A',adj=0,cex=1.5)
hist(ind_exp_filt$ad_sig, main="", ylab = "Adhesion Count", xlab = "Normalized Average Paxillin Intensity")
hist(ind_exp_filt$ax, main="", ylab = "Adhesion Count",  xlab = "Axial Ratio")

#plot(ind_exp_filt$cent_dist, ind_exp_filt$ad_sig, xlab = expression(paste('Adhesion Area (', symbol("m"), m^2, ')',sep='')), ylab = "Normalized Average Paxillin Intensity", pch=19, cex=0.25)
#mtext('B',adj=0,cex=1.5)
#plot(ind_exp_filt$cent_dist,ind_exp_filt$ax, xlab = expression(paste('Adhesion Area (', symbol("m"), m^2, ')',sep='')), ylab = "Axial Ratio", pch=19, cex=0.25)

smoothScatter(ind_exp_filt$cent_dist, ind_exp_filt$ad_sig, xlab = expression(paste('Distance from Cell Edge (', symbol("m"), 'm)',sep='')), ylab = "Normalized Average Paxillin Intensity", nrpoints = 0)
mtext('B',adj=0,cex=1.5)
smoothScatter(ind_exp_filt$cent_dist,ind_exp_filt$ax, xlab = expression(paste('Distance from Cell Edge (', symbol("m"), 'm)',sep='')), ylab = "Axial Ratio", nrpoints = 0)

mtext('C',adj=0,cex=1.5)
graphics.off()

########################################
#Comparing Area/Kinetics
########################################
pdf(file.path(out_folder,'area_vs_intensity_rsq.pdf'),height=12)
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

pdf(file.path(out_folder,'area_vs_intensity_slope.pdf'))
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

########################################
#Kinetics Figure
########################################
pdf(file.path(out_folder,'kinetics.pdf'),height=9)
layout(rbind(c(1,2),c(3,4),c(5)))

par(bty='n', mar=c(5,4.2,2,0))
plot_ad_seq(results[[1]],675);
mtext('A',adj=0,cex=1.5)
plot_ad_seq(results[[1]],675,type='disassembly');
mtext('B',adj=0,cex=1.5)

par(mar=c(4.5,2.5,2,0))
hist(results_nofilt$a$R, main='Assembly', ylab='', xlab=paste('Adjusted R Squared Values (n=',length(results_nofilt$a$R_sq),')', sep=''), freq=TRUE)
mtext('C',adj=0,cex=1.5)
hist(results_nofilt$dis$R,main='Disassembly',ylab='',xlab=paste('Adjusted R Squared Values (n=',length(results_nofilt$d$R_sq),')', sep=''), freq=TRUE)
mtext('D',adj=0,cex=1.5)

par(mar=c(2.6,5,2,0))
boxplot_with_points(list(results_filt$a$slope,results_filt$dis$slope), names=c(paste('Assembly (n=',length(results_filt$a$R_sq),')', sep=''), paste('Disassembly (n=',length(results_filt$dis$R_sq),')', sep='')), boxwex=0.6, ylab=expression(paste('Rate (',min^-1,')',sep='')))
mtext('E',adj=0,cex=1.5)

graphics.off()

########################################
#Spacial Figure
########################################
pdf(file.path(out_folder,'spacial.pdf'), height=9)
par(bty='n',mar=c(4.2,4.1,2,0.2))
layout(rbind(c(1,2),c(3,4),c(5)), width=c(1,1,0.5), height=c(1,1,2))
start_props = hist(results_filt$a$edge_dist[! is.na(results_filt$a$stable_lifetime)], xlab=expression(paste('Distance from Edge at Birth (', mu, 'm)', sep='')), main = '', ylab = '')
mtext('A',adj=0,cex=1.5)

hist(results_filt$dis$edge_dist[! is.na(results_filt$dis$stable_lifetime)] ,xlab=expression(paste('Distance from Edge at Death (', mu, 'm)', sep='')), main = '', ylab = '')
mtext('B',adj=0,cex=1.5)

plot(results_filt$a$edge_dist[! is.na(results_filt$a$stable_lifetime)], results_filt$a$stable_lifetime[! is.na(results_filt$a$stable_lifetime)], ylab='Stable Lifetime (min)',xlab=expression(paste('Distance from Edge at Birth (', mu, 'm)', sep='')), xlim=range(start_props$breaks))
mtext('C',adj=0,cex=1.5)

plot(results_filt$dis$edge_dist[! is.na(results_filt$dis$stable_lifetime)], results_filt$dis$stable_lifetime[! is.na(results_filt$dis$stable_lifetime)], ylab='Stable Lifetime (min)', xlab=expression(paste('Distance from Edge at Death (', mu, 'm)', sep='')), xlim=range(start_props$breaks))
mtext('D',adj=0,cex=1.5)

plot(results_nofilt$j$b, results_nofilt$j$d, xlab=expression(paste('Distance from Edge at Birth (', mu, 'm)', sep='')), ylab=expression(paste('Distance from Edge at Death (', mu, 'm)', sep='')), pch=20, cex=0.75)
birth_vs_death_model <- lm(death_dist ~ birth_dist, data=results_nofilt$joint)
abline(birth_vs_death_model, col='green', lwd = 2)
model_summary <- summary(birth_vs_death_model)
x_data <- data.frame(birth_dist = seq(min(results_nofilt$j$b),max(results_nofilt$j$b),by=0.01))
conf_int = predict(birth_vs_death_model, x_data, interval="confidence", level=0.95)
lines(x_data$birth_dist, conf_int[,2], col='red', lty=2, lwd = 2)
lines(x_data$birth_dist, conf_int[,3], col='red', lty=2, lwd = 2)
segments(0,0,max(results_nofilt$j$b),max(results_nofilt$j$b), col='blue', lty=4, lwd = 2)

x_eq_y_point = coef(birth_vs_death_model)[[1]]/(1-coef(birth_vs_death_model)[[2]])
#segments(x_eq_y_point,0,x_eq_y_point,x_eq_y_point, col='blue')
text(x_eq_y_point,25,paste('y=',sprintf('%.03f',coef(birth_vs_death_model)[[2]]),'x + ',sprintf('%.03f',coef(birth_vs_death_model)[[1]])), srt=atan(coef(birth_vs_death_model)[[2]])*(360/(2*pi)), cex=1.5)
mtext('E',adj=0,cex=1.5)
graphics.off()