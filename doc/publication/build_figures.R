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
#exp_dirs <- Sys.glob('../../results/lin_region_variation/S178A/no_log_trans/*/adhesion_props/models/')
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
area_nofilt = filter_results(area, primary_filter_results = results, min_R_sq = -Inf, max_p_val = Inf)

results_S_props = gather_general_props(results_S)
ind_exp_filt_S = gather_single_image_props(ind_results_S)
results_filt = filter_results(results)
results_S_filt = filter_results(results_S)
area_filt = filter_results(area, primary_filter_results = results)

#all_data = list(results_props = results_props, ind_exp_filt = ind_exp_filt, 
#	results_nofilt = results_nofilt, results_S_nofilt = results_S_nofilt, area_nofilt = area_nofilt,
#	results_S_props = results_S_props, ind_exp_filt_S = ind_exp_filt_S, results_filt = results_filt,
#	results_S_filt = results_S_filt, area_filt = area_filt)
#save(all_data, file='all_data.Rdata')

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
svg(file.path(out_folder,'statics','statics.svg'),width=12)
if (make_pdfs) {
	pdf(file.path(out_folder,'statics','statics.pdf'),width=12)
}
layout(rbind(c(1,2,3),c(4,5,6)),heights=c(1,0.5))
par(bty='n', mar=c(5,4.2,2,0.1))

plot.new()
mtext('A',adj=-.15,side=3,line=-0.5,cex=1.5)

plot.new()
mtext('B',adj=-.15,side=3,line=-0.5,cex=1.5)

plot.new()
mtext('C',adj=-.15,side=3,line=-0.5,cex=1.5)

area_data = ind_exp_filt$Area[ind_exp_filt$Area < 5];
hist_props = hist(area_data,plot=FALSE);

hist(area_data, main="", ylab = "Adhesion Count", xlab = expression(paste('Adhesion Area (', symbol("m"), m^2, ')',sep='')))
mtext('D',adj=-.15,side=3,line=-0.5,cex=1.5)
hist(ind_exp_filt$ad_sig, main="", ylab = "Adhesion Count", xlab = "Normalized Average Paxillin Intensity")
hist(ind_exp_filt$ax[ind_exp_filt$ax < 8], main="", ylab = "Adhesion Count",  xlab = "Axial Ratio")

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

par(bty='n', mar=c(4,4.2,1.1,0))

plot.new()
mtext('A',adj=-.17,side=3,line=-0.5,cex=1.5)

plot_ad_seq(results[[1]],675,type='overall')
mtext('B',adj=-.17,side=3,line=-0.5,cex=1.5)

plot_ad_seq(results[[1]],675, main = 'Assembly');
text(3,0.65,pos=3,expression(paste(R^2,' = 0.949')))
text(3,0.6,pos=3,adj=0,'Slope = 0.106')
mtext('C',adj=-.17,side=3,line=-0.5,cex=1.5)

plot_ad_seq(results[[1]],675,type='disassembly', main = 'Disassembly')
text(3,0.33,pos=3,expression(paste(R^2,' = 0.961')))
text(3,0.3,pos=3,adj=0,'Slope = 0.035')
mtext('D',adj=-.17,side=3,line=-0.5,cex=1.5)

par(bty='n', mar=c(2.1,4.2,1.1,0))

boxplot_with_points(list(results_filt$a$slope,results_filt$dis$slope), 
					names=c('Assembly', 'Disassembly'), boxwex=0.6, 
					ylab=expression(paste('Rate (',min^-1,')',sep='')))
#95% confidence intervals on the mean from Webb 2004
#segments(1.4,0.04,1.4,0.2,lwd=2)
#segments(1.35,0.12,1.45,0.12,lwd=2)
#segments(2.4,0.08,2.4,.088+0.004*2,lwd=2)
#segments(2.35,0.088,2.45,.088,lwd=2)
mtext('E',adj=-0.075,side=3,line=-0.5,cex=1.5)

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

start_props = hist(results_filt$a$edge_dist, 					   xlab=expression(paste('Distance from Edge at Birth (', mu, 'm) n=309', sep='')), 
				   main = '', ylab = '# of Focal Adhesions')
mtext('A',adj=-.2,side=3,line=-0.5,cex=1.5)

death_props = hist(results_filt$dis$edge_dist,
     xlab=expression(paste('Distance from Edge at Death (', mu, 'm) n=470', sep='')), 
     main = '', ylab = '# of Focal Adhesions')
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


########################################
#Forced Axis
########################################

svg(file.path(out_folder,'spacial','spacial_axes.svg'))
if (make_pdfs) {
	pdf(file.path(out_folder,'spacial','spacial_axes.pdf'))
}
par(bty='n',mar=c(4.2,4.1,2,0))
layout(rbind(c(1,2),c(3,4)))

start_props = hist(results_filt$a$edge_dist, 					   xlab=expression(paste('Distance from Edge at Birth (', mu, 'm) n=309', sep='')), 
				   main = '', ylab = '# of Focal Adhesions',			       axes=FALSE, xlim=c(0,30))
axis(2)
axis(1, at=seq(0,30, by=5))
mtext('A',adj=-.2,side=3,line=-0.5,cex=1.5)

death_props = hist(results_filt$dis$edge_dist,
     xlab=expression(paste('Distance from Edge at Death (', mu, 'm) n=470', sep='')), 
     main = '', ylab = '# of Focal Adhesions',
     axes=FALSE, xlim=c(0,30))
axis(2)
axis(1, at=seq(0,30, by=5))
mtext('B',adj=-.2,side=3,line=-0.5,cex=1.5)

par(bty='n',mar=c(4.2,4.1,2,0.4))

plot(results_filt$a$edge_dist, 
	 results_filt$a$slope, 
	 ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')),
	 xlab=expression(paste('Distance from Edge at Birth (', mu, 'm)', sep='')),
	 axes=FALSE, xlim=c(0,30))
axis(2)
axis(1, at=seq(0,30, by=5))
mtext('C',adj=-.2,side=3,line=-0.5,cex=1.5)

par(bty='n',mar=c(4.2,4.1,2,0.95))

plot(results_filt$d$edge_dist, 
	 results_filt$d$slope, 
	 ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')),
	 xlab=expression(paste('Distance from Edge at Death (', mu, 'm)', sep='')),
     axes=FALSE, xlim=c(0,30))
axis(2)
axis(1, at=seq(0,30, by=5))
mtext('D',adj=-.2,side=3,line=-0.5,cex=1.5)
graphics.off()

########################################
#Supplemental
########################################

svg(file.path(out_folder,'supplemental','spacial_nofilt.svg'))
if (make_pdfs) {
	pdf(file.path(out_folder,'supplemental','spacial_nofilt.svg'))
}
par(bty='n',mar=c(4.2,4.1,2,0))
layout(rbind(c(1,2),c(3,4)))

start_props = hist(results_nofilt$a$edge_dist, 					   xlab=expression(paste('Distance from Edge at Birth (', mu, 'm) n=309', sep='')), 
				   main = '', ylab = '# of Focal Adhesions')
mtext('A',adj=-.2,side=3,line=-0.5,cex=1.5)

death_props = hist(results_nofilt$dis$edge_dist,
     xlab=expression(paste('Distance from Edge at Death (', mu, 'm) n=470', sep='')), 
     main = '', ylab = '# of Focal Adhesions')
mtext('B',adj=-.2,side=3,line=-0.5,cex=1.5)

par(bty='n',mar=c(4.2,4.1,2,0.4))

plot(results_nofilt$a$edge_dist, 
	 results_nofilt$a$slope, 
	 ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')),
	 xlab=expression(paste('Distance from Edge at Birth (', mu, 'm)', sep='')))
mtext('C',adj=-.2,side=3,line=-0.5,cex=1.5)

par(bty='n',mar=c(4.2,4.1,2,0.95))

plot(results_nofilt$d$edge_dist, 
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
conf_int = predict(birth_vs_death_model, x_data, interval="confidence", level=0.95)
lines(x_data$birth_dist, conf_int[,2], col='red', lty=2, lwd = 3)
lines(x_data$birth_dist, conf_int[,3], col='red', lty=2, lwd = 3)

segments(0,0,max(results_nofilt$j$b),max(results_nofilt$j$b), col='blue', lty=4, lwd = 3)

graphics.off()
print('Done with Spacial')

############################################################
#Comparing S178A to Wild-type
############################################################

########################################
#Statics
########################################
wt = c()
for (i in 1:length(cell_area)) {
	wt  = c(wt, as.numeric(ad_area[[i]]/cell_area[[i]]))
}

s178a = c()
for (i in 1:length(cell_area_S)) {
	s178a  = c(s178a, as.numeric(ad_area_S[[i]]/cell_area_S[[i]]))
}

stage_data <- gather_stage_lengths(results_filt, results_S_filt)
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
par(bty='n',mar=c(2,4.1,1,0.2))
barplot(bar_lengths, names=c('Wild-type','S178A'), 
		ylab='Time (min)', width=matrix(0.3,3,2), xlim=c(0,1),
		legend=c('Assembly','Stability','Disassembly'),ylim = c(0,max(err_bars)))
errbar(err_bars[,1],err_bars[,2],err_bars[,3],err_bars[,4],add=TRUE,cex=0.0001)
graphics.off()

#95% confidence interval notation
#text(conf_ints[1,1], bar_lengths[1,1]/2,
#	 substitute( (x-y), list(x = sprintf('%.01f',conf_ints[1,3]), 
#	 						 y = sprintf('%.01f',conf_ints[1,4]))))
#text(conf_ints[2,1], bar_lengths[1,1] + bar_lengths[2,1]/2,
#	 substitute( (x-y), list(x = sprintf('%.01f',conf_ints[2,3]), 
#	 						 y = sprintf('%.01f',conf_ints[2,4]))))
#text(conf_ints[3,1], bar_lengths[1,1] + bar_lengths[2,1] + bar_lengths[3,1]/2,
#	 substitute( (x-y), list(x = sprintf('%.01f',conf_ints[3,3]), 
#	 						 y = sprintf('%.01f',conf_ints[3,4]))))
#	 						 
#text(conf_ints[4,1], bar_lengths[1,2]/2,
#	 substitute( (x-y), list(x = sprintf('%.01f',conf_ints[4,3]), 
#	 						 y = sprintf('%.01f',conf_ints[4,4]))))
#text(conf_ints[5,1], bar_lengths[1,2] + bar_lengths[2,2]/2,
#	 substitute( (x-y), list(x = sprintf('%.01f',conf_ints[5,3]), 
#	 						 y = sprintf('%.01f',conf_ints[5,4]))))
#text(conf_ints[6,1], bar_lengths[1,2] + bar_lengths[2,2] + bar_lengths[3,2]/2,
#	 substitute( (x-y), list(x = sprintf('%.01f',conf_ints[6,3]), 
#	 						 y = sprintf('%.01f',conf_ints[6,4]))))

#Plus/minus notation
#text(conf_ints[1,1], bar_lengths[1,1]/2,
#	 substitute( x %+-% y, list(x = sprintf('%.01f',bar_lengths[1,1]), 
#	 							  y = sprintf('%.01f',mean(abs(conf_ints[1,3:4] - bar_lengths[1,1]))))))
#text(conf_ints[2,1], bar_lengths[1,1] + bar_lengths[2,1]/2,
#	 substitute( x %+-% y, list(x = sprintf('%.01f',bar_lengths[2,1]), 
#	 							  y = sprintf('%.01f',mean(abs(conf_ints[2,3:4] - bar_lengths[2,1]))))))
#text(conf_ints[3,1], bar_lengths[1,1] + bar_lengths[2,1] + bar_lengths[3,1]/2,
#	 substitute( x %+-% y, list(x = sprintf('%.01f',bar_lengths[3,1]), 
#	 							  y = sprintf('%.01f',mean(abs(conf_ints[3,3:4] - bar_lengths[3,1]))))))
#
#text(conf_ints[4,1], bar_lengths[1,2]/2,
#	 substitute( x %+-% y, list(x = sprintf('%.01f',bar_lengths[1,2]), 
#	 							  y = sprintf('%.01f',mean(abs(conf_ints[4,3:4] - bar_lengths[1,2]))))))
#text(conf_ints[5,1], bar_lengths[1,2] + bar_lengths[2,2]/2,
#	 substitute( x %+-% y, list(x = sprintf('%.01f',bar_lengths[2,2]), 
#	 							  y = sprintf('%.01f',mean(abs(conf_ints[5,3:4] - bar_lengths[2,2]))))))
#text(conf_ints[6,1], bar_lengths[1,2] + bar_lengths[2,2] + bar_lengths[3,1]/2,
#	 substitute( x %+-% y, list(x = sprintf('%.01f',bar_lengths[3,2]), 
#	 							  y = sprintf('%.01f',mean(abs(conf_ints[6,3:4] - bar_lengths[3,2]))))))

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
				    names=c('Wild-type','S178A'),					ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')))
#segments(1,0.145,1,0.149)
#segments(1,0.149,2,0.149)
#segments(2,0.149,2,0.145)
#text(1.5,0.153,'***', cex=1.5)
mtext('A',adj=-.25,side=3,line=0.5,cex=1.5)	    

boxplot_with_points(list(results_filt$dis$slope,results_S_filt$dis$slope), 
				    names=c('Wild-type','S178A'), 
				    ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')))
#segments(1,0.1,1,0.11)
#segments(1,0.149,2,0.149)
#segments(2,0.149,2,0.145)
#text(1.5,0.153,'***', cex=1.5)
				    
mtext('B',adj=-.25,side=3,line=0.5,cex=1.5)	    

boxplot_with_points(list(results_filt$as$edge_dist,results_S_filt$as$edge_dist), 
					names=c('Wild-type','S178A'), 
					ylab=expression(paste('Distance from Edge at Birth (',mu,'m)',sep='')))
mtext('C',adj=-.25,side=3,line=0.5,cex=1.5)	    

boxplot_with_points(list(results_filt$dis$edge_dist,results_S_filt$dis$edge_dist), 
					names=c('Wild-type','S178A'), 
					ylab=expression(paste('Distance from Edge at Death (',mu,'m)',sep='')))
mtext('D',adj=-.25,side=3,line=0.5,cex=1.5)	    

graphics.off()

print('Done with S178A Comparisons')