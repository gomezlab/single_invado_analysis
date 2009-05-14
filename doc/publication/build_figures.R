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
svg(file.path(out_folder,'image_processing','image_processing.svg'),height=9)
#pdf(file.path(out_folder,'image_processing','image_processing.pdf'),height=9)
layout(rbind(c(1,1,1,2,2,2),c(3,3,4,4,5,5),c(6,6,6,7,7,7)),heights=c(1,0.75,1))
par(bty='n', mar=c(5,4.2,2,0.1))

plot.new()
mtext('A',adj=-.2,side=3,line=-0.5,cex=1.5)

plot.new()
mtext('B',adj=-.2,side=3,line=-0.5,cex=1.5)

hist(ind_exp_filt$Area, main="", ylab = "Adhesion Count", xlab = expression(paste('Adhesion Area (', symbol("m"), m^2, ')',sep='')))
mtext('C',adj=-.35,side=3,line=-0.5,cex=1.5)
hist(ind_exp_filt$ad_sig, main="", ylab = "Adhesion Count", xlab = "Normalized Average Paxillin Intensity")
hist(ind_exp_filt$ax, main="", ylab = "Adhesion Count",  xlab = "Axial Ratio")

plot.new()
mtext('D',adj=-.2,side=3,line=-0.5,cex=1.5)

plot_ad_seq(results[[1]],675,type='overall')
mtext('E',adj=-.2,side=3,line=-0.5,cex=1.5)

#smoothScatter(ind_exp_filt$cent_dist, ind_exp_filt$ad_sig, xlab = expression(paste('Distance from Cell Edge (', symbol("m"), 'm)',sep='')), ylab = "Normalized Average Paxillin Intensity", nrpoints = 0)
#mtext('B',adj=0,cex=1.5)
#smoothScatter(ind_exp_filt$cent_dist,ind_exp_filt$ax, xlab = expression(paste('Distance from Cell Edge (', symbol("m"), 'm)',sep='')), ylab = "Axial Ratio", nrpoints = 0)
#
#mtext('C',adj=0,cex=1.5)
graphics.off()
print('Done with Static Properties')

########################################
#Kinetics Figure
########################################
svg(file.path(out_folder,'kinetics','kinetics.svg'),height=9)
#pdf(file.path(out_folder,'kinetics','kinetics.pdf'),height=9)
layout(rbind(c(1,2),c(3,4),c(5,5)))

par(bty='n', mar=c(4,4.2,2,0))

plot_ad_seq(results[[1]],675, main = 'Assembly');
text(3,0.65,pos=3,expression(paste(R^2,' = 0.949')))
text(3,0.6,pos=3,adj=0,'Slope = 0.106')
mtext('A',adj=-.17,side=3,line=-0.5,cex=1.5)

plot_ad_seq(results[[1]],675,type='disassembly', main = 'Disassembly')
text(36,0.33,pos=3,expression(paste(R^2,' = 0.961')))
text(36,0.3,pos=3,adj=0,'Slope = 0.035')
mtext('B',adj=-.17,side=3,line=-0.5,cex=1.5)

hist(results_nofilt$a$R, main='Assembly' ,xlab=paste('Adjusted R Squared Values (n=',length(results_nofilt$a$R_sq),')', sep=''), freq=TRUE, ylab='# of Focal Adhesions')
mtext('C',adj=-.17,side=3,line=-0.5,cex=1.5)

hist(results_nofilt$dis$R,main='Disassembly', xlab=paste('Adjusted R Squared Values (n=',length(results_nofilt$d$R_sq),')', sep=''), freq=TRUE, ylab='# of Focal Adhesions')
mtext('D',adj=-.2,side=3,line=-0.5,cex=1.5)

par(bty='n', mar=c(2.1,4.2,2,0))

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

print('Done with Kinetics')

########################################
#Spacial Figure
########################################
svg(file.path(out_folder,'spacial','spacial.svg'), height=9)
pdf(file.path(out_folder,'spacial','spacial.pdf'), height=9)
par(bty='n',mar=c(4.2,4.1,2,0.2))
layout(rbind(c(1,2),c(3,4),c(5,6),c(7)), width=c(1,1,0.5), height=c(1,1,2))


start_props = hist(results_filt$a$edge_dist, 				   xlab=expression(paste('Distance from Edge at Birth (', mu, 'm)', sep='')), 
				   main = '', ylab = '')
mtext('A',adj=0,cex=1.5)


hist(results_filt$dis$edge_dist,
     xlab=expression(paste('Distance from Edge at Death (', mu, 'm)', sep='')), 
     main = '', ylab = '')
mtext('B',adj=0,cex=1.5)

plot(results_filt$a$edge_dist, 
	 results_filt$a$slope, 
	 ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')),
	 xlab=expression(paste('Distance from Edge at Birth (', mu, 'm)', sep='')))
mtext('C',adj=0,cex=1.5)


plot(results_filt$d$edge_dist, 
	 results_filt$d$slope, 
	 ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')),
	 xlab=expression(paste('Distance from Edge at Death (', mu, 'm)', sep='')))
mtext('D',adj=0,cex=1.5)

plot(results_nofilt$a$edge_dist[! is.na(results_nofilt$a$longevity)],
	 results_nofilt$a$longevity[! is.na(results_nofilt$a$longevity)], 
	 ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')),
	 xlab=expression(paste('Distance from Edge at Birth (', mu, 'm)', sep='')))
mtext('E',adj=0,cex=1.5)

plot(results_filt$d$edge_dist, 
	 results_filt$d$R_sq, 
	 ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')),
	 xlab=expression(paste('Distance from Edge at Birth (', mu, 'm)', sep='')))
mtext('E',adj=0,cex=1.5)


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

########################################
#Comparing S178A to Wild-type
########################################
svg(file.path(out_folder,'S178A','S178A_vs_wild-type.svg'))
#pdf(file.path(out_folder,'S178A','S178A_vs_wild-type.pdf'))
layout(rbind(c(1,2),c(3,4)))
par(bty='n', mar=c(2,4,2,0.1))

#boxplot_with_points(list(ind_exp_filt$Area,ind_exp_filt_S$Area), inc.points = FALSE,
#				    names=c('Wild-type','S178A'),
#				    ylab = expression(paste('Adhesion Area (', symbol("m"), m^2, ')',sep='')))
#mtext('A',adj=0,cex=1.5)

#boxplot_with_points(list(ind_exp_filt$axial,ind_exp_filt_S$axial), inc.points = FALSE,
#				    names=c('Wild-type','S178A'),
#				    ylab = 'Axial Ratio')
#mtext('B',adj=0,cex=1.5)

boxplot_with_points(list(results_filt$as$slope,results_S_filt$as$slope), 
				    names=c('Wild-type','S178A'),					ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')))
mtext('A',adj=-.25,side=3,line=-0.5,cex=1.5)	    

boxplot_with_points(list(results_filt$dis$slope,results_S_filt$dis$slope), 
				    names=c('Wild-type','S178A'), 
				    ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')))
mtext('B',adj=-.25,side=3,line=-0.5,cex=1.5)	    

boxplot_with_points(list(results_filt$as$edge_dist,results_S_filt$as$edge_dist), 
					names=c('Wild-type','S178A'), 
					ylab=expression(paste('Distance from Edge at Birth (',mu,'m)',sep='')))
mtext('C',adj=-.25,side=3,line=-0.5,cex=1.5)	    

boxplot_with_points(list(results_filt$dis$edge_dist,results_S_filt$dis$edge_dist), 
					names=c('Wild-type','S178A'), 
					ylab=expression(paste('Distance from Edge at Death (',mu,'m)',sep='')))
mtext('D',adj=-.25,side=3,line=-0.5,cex=1.5)	    

graphics.off()

########################################
#Comparing S178A to Wild-type
########################################
wt = c()
for (i in 1:length(cell_area)) {
	wt  = c(wt, as.numeric(ad_area[[i]]/cell_area[[i]]))
}

s178a = c()
for (i in 1:length(cell_area_S)) {
	s178a  = c(s178a, as.numeric(ad_area_S[[i]]/cell_area_S[[i]]))
}
