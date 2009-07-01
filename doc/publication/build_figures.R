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

exp_dirs <- Sys.glob('../../results/focal_adhesions/*/adhesion_props/models/')
#exp_dirs <- Sys.glob('../../results/lin_region_variation/FA/6/*/adhesion_props/models/')
#exp_dirs <- Sys.glob('../../results/lin_region_variation/FA/8/*/adhesion_props/models/')
#exp_dirs <- Sys.glob('../../results/lin_region_variation/FA/no_log_trans/*/adhesion_props/models/')
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
raw_data$results = load_results(exp_dirs,file.path('intensity_model.Rdata'))
raw_data$area = load_results(exp_dirs,file.path('area_model.Rdata'))
raw_data$ind_results <- load_data_files(exp_dirs, file.path('..','individual_adhesions.csv'), 						headers=T, debug=FALSE, inc_exp_names=FALSE);

single_props$fa$cell_area <- load_data_files(exp_dirs, 
	file.path('..','single_props','Cell_size.csv'), 
	headers=F, debug=FALSE, inc_exp_names=FALSE);
single_props$fa$ad_area <- load_data_files(exp_dirs,
	file.path('..','single_props','Ad_size.csv'), 
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
single_props$fa$outside_int <- load_data_files(exp_dirs,
	file.path('..','single_props','Outside_mean_intensity.csv'), 
	headers=F, debug=FALSE, inc_exp_names=FALSE);

exp_dirs_S <- Sys.glob('../../results/S178A/*/adhesion_props/models/')
#exp_dirs_S <- Sys.glob('../../results/lin_region_variation/S178A/6/*/adhesion_props/models/')
#exp_dirs_S <- Sys.glob('../../results/lin_region_variation/S178A/8/*/adhesion_props/models/')
#exp_dirs <- Sys.glob('../../results/lin_region_variation/S178A/no_log_trans/*/adhesion_props/models/')
exp_dirs_S <- exp_dirs_S[file_test('-d',exp_dirs_S)]
raw_data$results_S = load_results(exp_dirs_S,file.path('intensity_model.Rdata'))
raw_data$area_S = load_results(exp_dirs_S,file.path('area_model.Rdata'))
raw_data$ind_results_S <- load_data_files(exp_dirs_S,file.path('..','individual_adhesions.csv'), headers=T, debug=FALSE, inc_exp_names=FALSE);

cell_area_S <- load_data_files(exp_dirs_S,file.path('..','Cell_size.csv'), headers=F, debug=FALSE, inc_exp_names=FALSE);
ad_area_S <- load_data_files(exp_dirs_S,file.path('..','Ad_size.csv'), headers=F, debug=FALSE, inc_exp_names=FALSE);
mean_int_S <- load_data_files(exp_dirs_S,file.path('..','Mean_intensity.csv'), headers=F, debug=FALSE, inc_exp_names=FALSE);

print('Done Loading Data')

########################################
#Result filtering
########################################
results_props = gather_general_props(raw_data$results)
ind_exp_filt = gather_single_image_props(raw_data$ind_results)
results_nofilt = filter_results(raw_data$results, min_R_sq = -Inf, max_p_val = Inf)
results_S_nofilt = filter_results(raw_data$results_S, min_R_sq = -Inf, max_p_val = Inf)

results_onlysignif = filter_results(raw_data$results, min_R_sq = -Inf, max_p_val = 0.05)
results_S_onlysignif = filter_results(raw_data$results_S, min_R_sq = -Inf, max_p_val = 0.05)

area_nofilt = filter_results(raw_data$area, primary_filter_results = raw_data$results, min_R_sq = -Inf, max_p_val = 0.05)

results_S_props = gather_general_props(raw_data$results_S)
ind_exp_filt_S = gather_single_image_props(raw_data$ind_results_S)
results_filt = filter_results(raw_data$results)
results_all_filt = filter_results(raw_data$all_merge_results)
results_S_filt = filter_results(raw_data$results_S)
area_filt = filter_results(raw_data$area, primary_filter_results = raw_data$results)

rm(raw_data)
gc()

print('Done Filtering Data')

################################################################################
#Plotting
################################################################################
out_folder = '../../doc/publication/figures/'
make_pdfs = FALSE;

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
#General Properties
########################################
svg(file.path(out_folder,'statics','statics.svg'))
if (make_pdfs) {
	pdf(file.path(out_folder,'statics','statics.pdf'),width=8)
}
layout(rbind(c(1,1,2,2,3,3),c(4,4,4,5,5,5),c(6,6,6,7,7,7)),heights=c(1,0.85,0.85))
par(bty='n', mar=c(0,4,1.6,0))

plot.new()
mtext('A',adj=-.31,side=3,line=0,cex=1.5)

plot.new()
mtext('B',adj=-.31,side=3,line=0,cex=1.5)

plot.new()
mtext('C',adj=-.31,side=3,line=0,cex=1.5)

par(bty='n', mar=c(5,4.2,2,0.1))
area_data = ind_exp_filt$Area[ind_exp_filt$Area < 5];

hist(area_data, main="", ylab = "Adhesion Count", 
	 xlab = expression(paste('Adhesion Area (', symbol("m"), m^2, ')',sep='')))
mtext('D',adj=-.2,side=3,line=-0.25,cex=1.5)
hist(ind_exp_filt$ad_sig, main="", ylab = "Adhesion Count", xlab = "Normalized Average Paxillin Intensity")
hist(ind_exp_filt$ax[ind_exp_filt$ax < 8], main="", ylab = "Adhesion Count",  xlab = "Axial Ratio")
hist(ind_exp_filt$cent_dist, main="", ylab = "Adhesion Count",  
	 xlab = expression(paste("Distance from Edge (", symbol("m"), m, ')',sep='')))

graphics.off()
print('Done with Static Properties')

########################################
#Kinetics Figure
########################################
svg(file.path(out_folder,'kinetics','kinetics.svg'),height=9)
if (make_pdfs) {
	pdf(file.path(out_folder,'kinetics','kinetics.pdf'),height=9)
}
layout(rbind(c(1,2),c(3,4),c(5,5)))

exp_one_only = load_results(exp_dirs[[1]],'corrected_intensity_model.Rdata')

par(bty='n', mar=c(4,4.2,1.1,0))

plot.new()
mtext('A',adj=-.19,side=3,line=-0.5,cex=1.5)

ad_num = 675
plot_ad_seq(exp_one_only,ad_num,type='overall', 
	phase_lengths=c(exp_one_only$assembly$length[ad_num],exp_one_only$disassembly$length[ad_num]))
mtext('B',adj=-.19,side=3,line=-0.5,cex=1.5)

plot_ad_seq(exp_one_only,675, main = 'Assembly');
text(3,0.65,pos=3,expression(paste(R^2,' = 0.949')))
text(3,0.6,pos=3,adj=0,'Slope = 0.106')
mtext('C',adj=-.19,side=3,line=-0.5,cex=1.5)

plot_ad_seq(exp_one_only,675,type='disassembly', main = 'Disassembly')
text(3,0.33,pos=3,expression(paste(R^2,' = 0.961')))
text(3,0.3,pos=3,adj=0,'Slope = 0.035')
mtext('D',adj=-.19,side=3,line=-0.5,cex=1.5)

par(bty='n', mar=c(2.1,4.2,1.1,0))

boxplot_with_points(list(results_filt$a$slope,results_filt$dis$slope), 
					names=c('Assembly', 'Disassembly'), boxwex=0.6, 
					ylab=expression(paste('Rate (',min^-1,')',sep='')))
#95% confidence intervals on the mean from Webb 2004
#segments(1.4,0.04,1.4,0.2,lwd=2)
#segments(1.35,0.12,1.45,0.12,lwd=2)
#segments(2.4,0.08,2.4,.088+0.004*2,lwd=2)
#segments(2.35,0.088,2.45,.088,lwd=2)
mtext('E',adj=-0.085,side=3,line=-0.5,cex=1.5)

graphics.off()

svg(file.path(out_folder,'supplemental','R_squared.svg'),width=12)
if (make_pdfs) {
	pdf(file.path(out_folder,'supplemental','R_squared.pdf'),width=12)
}
layout(cbind(1,2))

par(bty='n', mar=c(4,4.2,2,0))

hist(results_nofilt$a$R, main='Assembly', freq=TRUE,
	 xlab=paste('Adjusted R Squared Values (n=',length(results_nofilt$a$R_sq),')', sep=''), 
	 ylab='# of Focal Adhesions')
mtext('A',adj=-.17,side=3,line=-0.5,cex=1.5)

hist(results_nofilt$dis$R,main='Disassembly', freq=TRUE,
	 xlab=paste('Adjusted R Squared Values (n=',length(results_nofilt$d$R_sq),')', sep=''), 
	 ylab='# of Focal Adhesions')
mtext('B',adj=-.2,side=3,line=-0.5,cex=1.5)
graphics.off()

print('Done with Kinetics')

########################################
#Spacial Figure
########################################
svg(file.path(out_folder,'spacial','spacial.svg'))
if (make_pdfs) {
	pdf(file.path(out_folder,'spacial','spacial.pdf'))
}
par(bty='n',mar=c(4.2,4.1,2,0))
layout(rbind(c(1,2),c(3,4)))

birth_props = hist(results_filt$a$edge_dist, 
				   xlab=expression(paste('Distance from Edge at Birth (', mu, 'm) n=309', sep='')), 
				   main = '', ylab = '# of Focal Adhesions', breaks=seq(0,26,by=2))
mtext('A',adj=-.2,side=3,line=-0.5,cex=1.5)

death_props = hist(results_filt$dis$edge_dist,
     xlab=expression(paste('Distance from Edge at Death (', mu, 'm) n=470', sep='')), 
     main = '', ylab = '# of Focal Adhesions', breaks=birth_props$breaks)
mtext('B',adj=-.2,side=3,line=-0.5,cex=1.5)

par(bty='n',mar=c(4.2,4.1,2,0.4))

plot(results_filt$a$edge_dist, 
	 results_filt$a$slope, 
	 ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')),
	 xlab=expression(paste('Distance from Edge at Birth (', mu, 'm)', sep='')))
mtext('C',adj=-.2,side=3,line=-0.5,cex=1.5)

par(bty='n',mar=c(4.2,4.1,2,0.95))

plot(results_filt$d$edge_dist, 
	 results_filt$d$slope, 
	 ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')),
	 xlab=expression(paste('Distance from Edge at Death (', mu, 'm)', sep='')))
mtext('D',adj=-.2,side=3,line=-0.5,cex=1.5)
graphics.off()

####################
#Supplemental
####################
svg(file.path(out_folder,'supplemental','spacial_nofilt.svg'))
if (make_pdfs) {
	pdf(file.path(out_folder,'supplemental','spacial_nofilt.pdf'))
}
par(bty='n',mar=c(4.2,4.1,2,0))
layout(rbind(c(1,2),c(3,4)))

start_props = hist(results_nofilt$a$edge_dist, 
				   xlab=expression(paste('Distance from Edge at Birth (', mu, 'm) n=1169', sep='')), 
				   main = '', ylab = '# of Focal Adhesions',breaks=seq(0,38,by=1))
mtext('A',adj=-.2,side=3,line=-0.5,cex=1.5)

death_props = hist(results_nofilt$dis$edge_dist,
     xlab=expression(paste('Distance from Edge at Death (', mu, 'm) n=1460', sep='')), 
     main = '', ylab = '# of Focal Adhesions',breaks=seq(0,38,by=1))
mtext('B',adj=-.2,side=3,line=-0.5,cex=1.5)

par(bty='n',mar=c(4.2,4.1,2,0.4))

smoothScatter(results_nofilt$a$edge_dist, 
	 results_nofilt$a$slope, 
	 ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')),
	 xlab=expression(paste('Distance from Edge at Birth (', mu, 'm)', sep='')))
mtext('C',adj=-.2,side=3,line=-0.5,cex=1.5)

par(bty='n',mar=c(4.2,4.1,2,0.95))

smoothScatter(results_nofilt$d$edge_dist, 
	 results_nofilt$d$slope, 
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
x_data <- data.frame(birth_dist = seq(min(results_nofilt$j$b),max(results_nofilt$j$b),by=0.01))
line_conf = predict(birth_vs_death_model, x_data, interval="confidence", level=0.95)
lines(x_data$birth_dist, line_conf[,2], col='red', lty=2, lwd = 3)
lines(x_data$birth_dist, line_conf[,3], col='red', lty=2, lwd = 3)

segments(0,0,max(results_nofilt$j$b),max(results_nofilt$j$b), col='blue', lty=4, lwd = 3)

graphics.off()
print('Done with Spacial')

############################################################
#Comparing S178A to Wild-type
############################################################

########################################
#Statics
########################################
cors = c()
for (i in 1:length(cell_area)) {	
	cors = c(cors, cor(as.numeric(ad_area[[i]]/cell_area[[i]]), as.numeric(mean_int[[i]])))
}
cors_S = c()
for (i in 1:length(cell_area_S)) {	
	cors_S = c(cors_S, cor(as.numeric(ad_area_S[[i]]/cell_area_S[[i]]), as.numeric(mean_int_S[[i]])))
}

stage_data <- gather_stage_lengths(results_onlysignif, results_S_onlysignif, debug=TRUE)
#stage_data_filt <- gather_stage_lengths(results_filt, results_S_filt, debug=TRUE)
bar_lengths = stage_data$bar_lengths
conf_ints = stage_data$conf_ints
conf_ints[1:3,1] = 0.22
conf_ints[4:6,1] = 0.57

err_bars = conf_ints
err_bars[1,3:4] = err_bars[1,3:4]
err_bars[2,3:4] = err_bars[2,3:4] + sum(bar_lengths[1,1])
err_bars[3,3:4] = err_bars[3,3:4] + sum(bar_lengths[1:2,1])
err_bars[4,3:4] = err_bars[4,3:4]
err_bars[5,3:4] = err_bars[5,3:4] + sum(bar_lengths[1,2])
err_bars[6,3:4] = err_bars[6,3:4] + sum(bar_lengths[1:2,2])

pdf(file.path(out_folder,'lifetimes','adhesion_phase_lifetimes.pdf'))
par(bty='n',mar=c(2,4,0,0))
barplot(bar_lengths, names=c('Wild-type','S178A'), 
		ylab='Time (min)', width=matrix(0.3,3,2), xlim=c(0,1),
		legend=c('Assembly','Stability','Disassembly'),ylim = c(0,max(err_bars)+3.5))
errbar(err_bars[,1],err_bars[,2],err_bars[,3],err_bars[,4],add=TRUE,cex=0.0001)

#Signifcance Bars
bar_length = 1;
sep_from_data = 0.5;

upper_left = c(err_bars[1,1], min(err_bars[1,3],err_bars[4,3]) - sep_from_data);
lower_right = c(err_bars[4,1], min(err_bars[1,3],err_bars[4,3]) - (sep_from_data + bar_length));
lines(c(upper_left[1],upper_left[1],lower_right[1], lower_right[1]),
	  c(upper_left[2],lower_right[2],lower_right[2],upper_left[2]))  
text(mean(c(upper_left[1],lower_right[1]))-0.005,lower_right[2]-sep_from_data,"***",cex=1.5)

upper_left = c(err_bars[3,1], max(err_bars[3,4],err_bars[6,4]) + sep_from_data + bar_length);
lower_right = c(err_bars[6,1], max(err_bars[6,4],err_bars[6,4]) + sep_from_data);
lines(c(upper_left[1],upper_left[1],lower_right[1], lower_right[1]),
	  c(lower_right[2],upper_left[2],upper_left[2],lower_right[2]))  
text(mean(c(upper_left[1],lower_right[1]))-0.005,upper_left[2]+sep_from_data,"*",cex=1.5)

graphics.off()

########################################
#Dynamics
########################################
svg(file.path(out_folder,'S178A','S178A_vs_wild-type.svg'))
if (make_pdfs) {
	pdf(file.path(out_folder,'S178A','S178A_vs_wild-type.pdf'))
}
layout(rbind(c(1,2),c(3,4)))
par(bty='n', mar=c(2,4,2.3,0.1))

boxplot_with_points(list(results_filt$as$slope,results_S_filt$as$slope), 
				    names=c('Wild-type','S178A'),					ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')),ylim=c(0,0.15))
mtext('A',adj=-.25,side=3,line=0.5,cex=1.5)	    

boxplot_with_points(list(results_filt$dis$slope,results_S_filt$dis$slope), 
				    names=c('Wild-type','S178A'), 
				    ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')), ylim=c(0,0.10))
mtext('B',adj=-.25,side=3,line=0.5,cex=1.5)

boxplot_with_points(list(results_filt$as$edge_dist,results_S_filt$as$edge_dist), 
					names=c('Wild-type','S178A'), ylim=c(0,40), 
					ylab=expression(paste('Distance from Edge at Birth (',mu,'m)',sep='')))
mtext('C',adj=-.25,side=3,line=0.5,cex=1.5)

boxplot_with_points(list(results_filt$dis$edge_dist,results_S_filt$dis$edge_dist), 
					names=c('Wild-type','S178A'), ylim=c(0,40),
					ylab=expression(paste('Distance from Edge at Death (',mu,'m)',sep='')))
mtext('D',adj=-.25,side=3,line=0.5,cex=1.5)	    

graphics.off()

###############
#Conf Levels
###############
require(boot)
boot_one = boot(results_filt$as$slope, function(data,indexes) mean(data[indexes],na.rm=T), 50000);
boot_two = boot(results_S_filt$as$slope, function(data,indexes) mean(data[indexes],na.rm=T), 50000);
find_p_val_from_bootstrap(boot_one, boot_two)
boot_two$t0/boot_one$t0

boot_one = boot(results_filt$dis$slope, function(data,indexes) mean(data[indexes],na.rm=T), 50000);
boot_two = boot(results_S_filt$dis$slope, function(data,indexes) mean(data[indexes],na.rm=T), 50000);
find_p_val_from_bootstrap(boot_one, boot_two)
boot_two$t0/boot_one$t0

boot_one = boot(results_filt$as$edge_dist, function(data,indexes) median(data[indexes],na.rm=T), 50000);
boot_two = boot(results_S_filt$as$edge_dist, function(data,indexes) median(data[indexes],na.rm=T), 50000);
find_p_val_from_bootstrap(boot_one, boot_two)
boot_two$t0/boot_one$t0

boot_one = boot(results_filt$dis$edge_dist, function(data,indexes) median(data[indexes],na.rm=T), 50000);
boot_two = boot(results_S_filt$dis$edge_dist, function(data,indexes) median(data[indexes],na.rm=T), 50000);
find_p_val_from_bootstrap(boot_one, boot_two)

print('Done with S178A Comparisons')