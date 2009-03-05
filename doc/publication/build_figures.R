rm(list = ls())
source('linear_regions.R')
library(lattice)
library(geneplotter)
library(Hmisc)

################################################################################
#Result loading
################################################################################
exp_dirs <- Sys.glob('../../results/focal_adhesions/*/adhesion_props/lin_time_series/')
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
results = load_results(exp_dirs,file.path('..','intensity_model.Rdata'))
corr_results = load_results(exp_dirs,file.path('..','corrected_intensity_model.Rdata'))
area = load_results(exp_dirs,file.path('..','area_model.Rdata'))
box_results = load_results(exp_dirs,file.path('..','box_model.Rdata'))

ind_results <- load_data_files(exp_dirs,file.path('..','individual_adhesions.csv'), T);

########################################
#Result filtering
########################################
general_props = gather_general_props(results)

results_nofilt = filter_mixed_results(results, corr_results, min_R_sq = -Inf, max_p_val = Inf)
corr_nofilt = filter_mixed_area(corr_results, corr_results, min_R_sq = -Inf, max_p_val = Inf)
area_nofilt = filter_mixed_area(area, corr_results, min_R_sq = -Inf, max_p_val = Inf)
box_nofilt = filter_mixed_area(box_results, box_results, min_R_sq = -Inf, max_p_val = Inf)

results_filt = filter_mixed_results(results, corr_results)
corr_filt = filter_mixed_area(corr_results, corr_results)
area_filt = filter_mixed_area(area, corr_results)
box_filt = filter_mixed_area(box_results, box_results)

################################################################################
#Plotting
################################################################################
out_folder = '../../doc/publication/figures/'

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
plot_ind_data = list()

for (i in 1:length(ind_results)) {
	res = ind_results[[i]]
	filt_by_area = res$Area >= min(res$Area)*3 & res$I_num == 1
	plot_ind_data$Area = c(plot_ind_data$Area, res$Area[filt_by_area]);
	plot_ind_data$ad_sig = c(plot_ind_data$ad_sig, res$Average_adhesion_signal[filt_by_area]);
	plot_ind_data$axial_r = c(plot_ind_data$axial_r, res$MajorAxisLength[filt_by_area]/res$MinorAxisLength[filt_by_area]);
	
	plot_ind_data$cent_dist = c(plot_ind_data$cent_dist, res$Centroid_dist_from_edge[filt_by_area]);
}

pdf(file.path(out_folder,'general_props.pdf'))
layout(cbind(c(1,1,2,2,3,3), c(0,4,4,5,5,0)))
par(bty='n', mar=c(5,4.2,2,0.1))
hist(plot_ind_data$Area, main="", ylab = "Adhesion Count", xlab = expression(paste('Adhesion Area (', symbol("m"), m^2, ')',sep='')))
mtext('A',adj=0,cex=1.5)
hist(plot_ind_data$ad_sig, main="", ylab = "Adhesion Count", xlab = "Normalized Average Paxillin Intensity")
hist(plot_ind_data$ax, main="", ylab = "Adhesion Count",  xlab = "Axial Ratio")

#plot(plot_ind_data$cent_dist, plot_ind_data$ad_sig, xlab = expression(paste('Adhesion Area (', symbol("m"), m^2, ')',sep='')), ylab = "Normalized Average Paxillin Intensity", pch=19, cex=0.25)
#mtext('B',adj=0,cex=1.5)
#plot(plot_ind_data$cent_dist,plot_ind_data$ax, xlab = expression(paste('Adhesion Area (', symbol("m"), m^2, ')',sep='')), ylab = "Axial Ratio", pch=19, cex=0.25)

smoothScatter(plot_ind_data$cent_dist, plot_ind_data$ad_sig, xlab = expression(paste('Distance from Cell Edge (', symbol("m"), 'm)',sep='')), ylab = "Normalized Average Paxillin Intensity", nrpoints = 0)
mtext('B',adj=0,cex=1.5)
smoothScatter(plot_ind_data$cent_dist,plot_ind_data$ax, xlab = expression(paste('Distance from Cell Edge (', symbol("m"), 'm)',sep='')), ylab = "Axial Ratio", nrpoints = 0)

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
plot_ad_seq(corr_results[[1]],675);
mtext('A',adj=0,cex=1.5)
plot_ad_seq(corr_results[[1]],675,type='disassembly');
mtext('B',adj=0,cex=1.5)

par(mar=c(4.5,2.5,2,0))
hist(results_nofilt$a$R, main='Assembly', ylab='', xlab=paste('Adjusted R Squared Values (n=',length(results_nofilt$a$R_sq),')', sep=''), freq=TRUE)
mtext('C',adj=0,cex=1.5)
hist(results_nofilt$dis$R,main='Disassembly',ylab='',xlab=paste('Adjusted R Squared Values (n=',length(results_nofilt$d$R_sq),')', sep=''), freq=TRUE)
mtext('D',adj=0,cex=1.5)

par(mar=c(2.6,5,2,0))
boxplot_with_points(list(results_filt$a$slope,results_filt$dis$slope), names=c(paste('Assembly (n=',length(results_filt$a$R_sq),')', sep=''), paste('Disassembly (n=',length(results_filt$dis$R_sq),')', sep='')), boxwex=0.6, ylab=expression(paste('Rate (',min^-1,')',sep='')))
mtext('E',adj=0,cex=1.5)

#par(bty='n', mar=c(5,4.2,2,0))
#plot_ad_seq(corr_results[[1]],672);
#mtext('A',adj=0,cex=1.5)
#plot_ad_seq(corr_results[[1]],672,type='disassembly');
#mtext('B',adj=0,cex=1.5)
#
#par(mar=c(4.5,2.5,2,0))
#hist(box_nofilt$a$R, main='Assembly', ylab='', xlab=paste('Adjusted R Squared Values (n=',length(box_nofilt$a$R_sq),')', sep=''), freq=TRUE)
#mtext('C',adj=0,cex=1.5)
#hist(box_nofilt$dis$R,main='Disassembly',ylab='',xlab=paste('Adjusted R Squared Values (n=',length(box_nofilt$d$R_sq),')', sep=''), freq=TRUE)
#mtext('D',adj=0,cex=1.5)
#
#par(mar=c(2.6,5,2,0))
#boxplot_with_points(list(box_filt$a$slope,box_filt$dis$slope), names=c(paste('Assembly (n=',length(box_filt$a$R_sq),')', sep=''), paste('Disassembly (n=',length(box_filt$dis$R_sq),')', sep='')), boxwex=0.6, ylab=expression(paste('Rate (',min^-1,')',sep='')))
#mtext('E',adj=0,cex=1.5)

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

########################################
#Edge Velocity
########################################

library(boot)

files_to_load = c('../single_lin.csv', '../edge_velo_pre_birth.csv', 
				  '../edge_velo_post_birth.csv', '../edge_velo_null.csv', 
				  '../edge_velo_post_birth.csv', '../edge_velo_pre_death.csv')

raw <- load_data_files(exp_dirs, files_to_load, c(T,F,F,F, F, F), debug=T);
raw <- lapply(raw, function(x) x[2:10])

lineage_data <- raw[[1]]
pre_birth_raw <- raw[[2]]
post_death_raw <- raw[[3]]
null_raw <- raw[[4]]
post_birth_raw <- raw[[5]]
pre_death_raw <- raw[[6]]

pdf(file.path(out_folder,'ring_plots.pdf'),height=9, width=9)
par(mfcol=c(3,3),bty='n',oma=c(0,0,2,0), mar=c(4.2,4.1,2,0.2))

for (k in 1:10) {
	pre_birth_rings = list()
	pre_birth_ring_values = c(k,k)

	for (i in 1:length(pre_birth_raw)) {
		pre_birth_rings[[i]] = list()
		for (j in 1:length(pre_birth_ring_values)) {
			this_pre_birth = pre_birth_raw[[i]];
		
			birth_filt = ! lineage_data[[i]]$split_birth_status &
						 ! is.na(this_pre_birth[,length(this_pre_birth)]) &
						 lineage_data[[i]]$longevity > 0
			if (j == 1) {
				birth_filt = birth_filt & lineage_data[[i]]$starting_edge_dist < pre_birth_ring_values[j]
			}
			if (j == length(pre_birth_ring_values)) {
				birth_filt = birth_filt & lineage_data[[i]]$starting_edge_dist >= pre_birth_ring_values[j]
			} 
			if (j != 1 & j != length(pre_birth_ring_values)) {
				birth_filt = birth_filt & 
							 lineage_data[[i]]$starting_edge_dist > pre_birth_ring_values[j] &
							 lineage_data[[i]]$starting_edge_dist <= pre_birth_ring_values[j + 1]
			}
	
			this_pre_birth = this_pre_birth[birth_filt,]
			
			pre_birth_rings[[i]]$means[[j]] = as.numeric(mean(this_pre_birth, na.rm=T))
	
			temp = find_col_conf_ints(this_pre_birth, boot.samp=100)
			pre_birth_rings[[i]]$upper[[j]] = temp$u
			pre_birth_rings[[i]]$lower[[j]] = temp$l
		}
	}

	cols = rainbow(length(pre_birth_rings$means))
	cols = c("red", "green")
	for (i in 1:length(pre_birth_rings)) {
		ylims = c(0,0)
		ylims[1] = min(unlist(pre_birth_rings[[i]]$lower), na.rm=T)
		ylims[2] = max(unlist(pre_birth_rings[[i]]$upper), na.rm=T)	
		for (j in 1:length(pre_birth_rings[[i]]$means)) {
		
			add_status = F;
			if (j > 1) {
				add_status = T;
			}
		
			errbar(1:length(pre_birth_rings[[i]]$means[[j]]), 				   pre_birth_rings[[i]]$means[[j]], 
		     	   pre_birth_rings[[i]]$upper[[j]], 
			       pre_birth_rings[[i]]$lower[[j]], 
			       col = cols[j], add = add_status, ylab='',xlab='', ylim=ylims)
			if (j == 1) {
				mtext(paste('Red < ',k,'; Green >= ',k, sep=''), outer =T)
			}
		}
	}
}
graphics.off()



pre_birth = list()
post_birth = list()
pre_death = list()
post_death = list()
null = list()
for (i in 1:length(pre_birth_raw)) {
	this_pre_birth = pre_birth_raw[[i]];
	this_post_death = post_death_raw[[i]];
	this_null_raw = null_raw[[i]];
	this_post_birth = post_birth_raw[[i]];
	this_pre_death = pre_death_raw[[i]];
	
	birth_filt = ! lineage_data[[i]]$split_birth_status & 
					 !is.na(this_pre_birth[,length(this_pre_birth)]) & 
					 lineage_data[[i]]$starting_edge_dist < Inf
	
	death_filt = as.logical(lineage_data[[i]]$death_status) 
	
	this_pre_birth = this_pre_birth[birth_filt,]
	pre_birth$data[[i]] = this_pre_birth
	this_post_birth = this_post_birth[birth_filt,]
	post_birth$data[[i]] = this_post_birth
	this_pre_death = this_pre_death[death_filt,] 
	pre_death$data[[i]] = this_pre_death
	this_post_death = this_post_death[death_filt,] 
	post_death$data[[i]] = this_post_death

	pre_birth$means[[i]] = as.numeric(mean(this_pre_birth, na.rm=T))
	post_birth$means[[i]] = as.numeric(mean(this_post_birth, na.rm=T))
	pre_death$means[[i]] = as.numeric(mean(this_pre_death, na.rm=T))
	post_death$means[[i]] = as.numeric(mean(this_post_death, na.rm=T))
	
	temp = find_col_conf_ints(this_pre_birth)
	pre_birth$upper[[i]] = temp$u
	pre_birth$lower[[i]] = temp$l

	temp = find_col_conf_ints(this_post_birth)
	post_birth$upper[[i]] = temp$u
	post_birth$lower[[i]] = temp$l

	temp = find_col_conf_ints(this_pre_death)
	pre_death$upper[[i]] = temp$u
	pre_death$lower[[i]] = temp$l

	temp = find_col_conf_ints(this_post_death)
	post_death$upper[[i]] = temp$u
	post_death$lower[[i]] = temp$l

	pre_birth$unlist[[i]] = unlist(this_pre_birth)
	pre_birth$unlist[[i]] = pre_birth$unlist[[i]][!is.na(pre_birth$unlist[[i]])]
	temp = t.test(pre_birth$unlist[[i]])
	pre_birth$all_l[[i]] = temp$conf.int[1]
	pre_birth$all_u[[i]] = temp$conf.int[2]

	post_birth$unlist[[i]] = unlist(this_post_birth)
	post_birth$unlist[[i]] = post_birth$unlist[[i]][!is.na(post_birth$unlist[[i]])]
	temp = t.test(post_birth$unlist[[i]])
	post_birth$all_l[[i]] = temp$conf.int[1]
	post_birth$all_u[[i]] = temp$conf.int[2]

	pre_death$unlist[[i]] = unlist(this_pre_death)
	pre_death$unlist[[i]] = pre_death$unlist[[i]][!is.na(pre_death$unlist[[i]])]
	temp = t.test(pre_death$unlist[[i]])
	pre_death$all_l[[i]] = temp$conf.int[1]
	pre_death$all_u[[i]] = temp$conf.int[2]
	
	post_death$unlist[[i]] = unlist(this_post_death)
	post_death$unlist[[i]] = post_death$unlist[[i]][!is.na(post_death$unlist[[i]])]
	temp = t.test(post_death$unlist[[i]])
	post_death$all_l[[i]] = temp$conf.int[1]
	post_death$all_u[[i]] = temp$conf.int[2]
	
	null$unlist[[i]] = unlist(this_null_raw)
	null$unlist[[i]] = null$unlist[[i]][!is.na(null$unlist[[i]])]
	temp = t.test(null$unlist[[i]])
	null$all_l[[i]] = temp$conf.int[1]
	null$all_u[[i]] = temp$conf.int[2]
}

pdf(file.path(out_folder,'edge_experiment_plots.pdf'))
par(mfcol=c(3,3),bty='n',mar=c(4.2,4.1,2,0.2), oma=c(0,0,2,0))
for (i in 1:length(pre_birth$means)) {
	errbar(1:length(pre_birth$means[[i]]), pre_birth$means[[i]], pre_birth$upper[[i]], pre_birth$lower[[i]], ylab = i, xlab = "")

#	segments(1, null$all_u[[i]], length(pre_birth$means[[i]]), null$all_u[[i]], col='red')
#	segments(1, null$all_l[[i]], length(pre_birth$means[[i]]), null$all_l[[i]], col='red')

	segments(1, pre_birth$all_u[[i]], length(pre_birth$means[[i]]), pre_birth$all_u[[i]], col='green')
	segments(1, pre_birth$all_l[[i]], length(pre_birth$means[[i]]), pre_birth$all_l[[i]], col='green')
	
	if (i == 1) {	
		mtext('Pre-birth', outer=T)
	}
}

for (i in 1:length(post_birth$means)) {
	errbar(1:length(post_birth$means[[i]]), post_birth$means[[i]], post_birth$upper[[i]], post_birth$lower[[i]], ylab = i, xlab = "")
	
#	segments(1, null$all_u[[i]], length(post_birth$means[[i]]), null$all_u[[i]], col='red')
#	segments(1, null$all_l[[i]], length(post_birth$means[[i]]), null$all_l[[i]], col='red')

	segments(1, post_birth$all_u[[i]], length(post_birth$means[[i]]), post_birth$all_u[[i]], col='green')
	segments(1, post_birth$all_l[[i]], length(post_birth$means[[i]]), post_birth$all_l[[i]], col='green')
	if (i == 1) {
		mtext('Post-birth', outer=T)
	}
}

for (i in 1:length(pre_death$means)) {
	errbar(1:length(pre_death$means[[i]]), pre_death$means[[i]], pre_death$upper[[i]], pre_death$lower[[i]], ylab = i, xlab = "")
	
#	segments(1, null$all_u[[i]], length(pre_death$means[[i]]), null$all_u[[i]], col='red')
#	segments(1, null$all_l[[i]], length(pre_death$means[[i]]), null$all_l[[i]], col='red')

	segments(1, pre_death$all_u[[i]], length(pre_birth$means[[i]]), pre_death$all_u[[i]], col='green')
	segments(1, pre_death$all_l[[i]], length(pre_birth$means[[i]]), pre_death$all_l[[i]], col='green')

	if (i == 1) {
		mtext('Pre-death', outer=T)
	}
}

for (i in 1:length(post_death$means)) {
	errbar(1:length(post_death$means[[i]]), post_death$means[[i]], post_death$upper[[i]], post_death$lower[[i]], ylab = i, xlab = "")

#	segments(1, null$all_u[[i]], length(post_death$means[[i]]), null$all_u[[i]], col='red')
#	segments(1, null$all_l[[i]], length(post_death$means[[i]]), null$all_l[[i]], col='red')

	segments(1, post_death$all_u[[i]], length(pre_birth$means[[i]]), post_death$all_u[[i]], col='green')
	segments(1, post_death$all_l[[i]], length(pre_birth$means[[i]]), post_death$all_l[[i]], col='green')
	
	if (i == 1) {
		mtext('Post-death', outer=T)
	}
}

#for (i in 1:length(pre_birth$means)) {
#	errbar(1:10, 
#		pre_birth$means[[i]][(length(pre_birth$means[[i]]) - 9):length(pre_birth$means[[i]])], 
#		pre_birth$upper[[i]][(length(pre_birth$means[[i]]) - 9):length(pre_birth$means[[i]])], 
#		pre_birth$lower[[i]][(length(pre_birth$means[[i]]) - 9):length(pre_birth$means[[i]])],
#		ylab = i, xlab = "")
#	segments(1, null$all_u[[i]], length(pre_birth$means[[i]]), null$all_u[[i]], col='red')
#	segments(1, null$all_l[[i]], length(pre_birth$means[[i]]), null$all_l[[i]], col='red')
#}
#
#for (i in 1:length(post_death$means)) {
#	errbar(1:10, post_death$means[[i]][1:10], post_death$upper[[i]][1:10], post_death$lower[[i]][1:10], ylab = i, xlab = "")
#	segments(1, null$all_u[[i]], length(post_death$means[[i]]), null$all_u[[i]], col='red')
#	segments(1, null$all_l[[i]], length(post_death$means[[i]]), null$all_l[[i]], col='red')
#}

for (i in 1:length(pre_birth$unlist)) {
	hist(pre_birth$unlist[[i]], main=i,xlab='Projected Edge Velocity')
}

graphics.off()

########################################
#Split Birth Filter ROC Curve
########################################

#sb_thresh = list()
#sb_thresh$level = c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.01)
#sb_thresh$no_filt_num = c(2001, 1967, 1885, 1814, 1739, 1642, 1591, 1535, 1496, 1437, 1377)
#sb_thresh$no_filt_neg = c(519, 492, 442, 395, 358, 307, 282, 252, 234, 201, 171)
#sb_thresh$no_filt_pos = sb_thresh$no_filt_num - sb_thresh$no_filt_neg
#sb_thresh$filt_num = c(378, 375, 368, 358, 350, 341, 336, 331, 329, 323, 317)
#sb_thresh$filt_neg = c(39, 36, 30, 23, 22, 17, 15, 12, 12, 8, 5)
#sb_thresh$filt_pos = sb_thresh$filt_num - sb_thresh$filt_neg
#
#no_filt_init_num = 2202;
#no_filt_init_neg = 611;
#filt_init_num = 395;
#filt_init_neg = 44;
#
#true_pos = no_filt_init_neg - sb_thresh$no_filt_neg
#false_pos = no_filt_init_num - no_filt_init_neg - sb_thresh$no_filt_pos
#false_neg = sb_thresh$no_filt_neg
#true_neg = sb_thresh$no_filt_pos
#
#sensitivity = true_pos/(true_pos + false_neg)
#specificity = true_neg/(true_neg + false_pos)
#
#plot(1 - specificity, sensitivity, xlim=c(0,1), ylim=c(0,1), typ='l')
#segments(0,0,1,1)
#
#true_pos_filt = filt_init_neg - sb_thresh$filt_neg
#false_pos_filt = filt_init_num - filt_init_neg - sb_thresh$filt_pos
#false_neg_filt = sb_thresh$filt_neg
#true_neg_filt = sb_thresh$filt_pos
#
#sensitivity = true_pos_filt/(true_pos_filt + false_neg_filt)
#specificity = true_neg_filt/(true_neg_filt + false_pos_filt)
#
#plot(1 - specificity, sensitivity, xlim=c(0,1), ylim=c(0,1))
#segments(0,0,1,1)
