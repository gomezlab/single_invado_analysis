source('linear_regions.R')

################################################################################
#Result loading
################################################################################
results = list()
corr_results = list()
shrunk_results = list()

exp_dirs <- Sys.glob('../../results/focal_adhesions/*/adhesion_props/lin_time_series/')
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]
results = load_results(exp_dirs,file.path('..','intensity_model.Rdata'))
corr_results = load_results(exp_dirs,file.path('..','corrected_intensity_model.Rdata'))
shrunk_results = load_results(exp_dirs,file.path('..','shrunk_intensity_model.Rdata'))

################################################################################
#Result filtering
################################################################################
results_nofilt = filter_results(results, needed_R_sq = -Inf)
corr_nofilt = filter_results(corr_results, needed_R_sq = -Inf)
shrunk_nofilt = filter_results(shrunk_results, needed_R_sq = -Inf)

results_filt = filter_results(results)
corr_filt = filter_results(corr_results)
shrunk_filt = filter_results(shrunk_results)

################################################################################
#Plotting
################################################################################
out_folder = '../../doc/publication/figures/'

#Dyanamics Figure
pdf(file.path(out_folder,'dynamics.pdf'),height=9)
layout(rbind(c(1,2),c(3,4),c(5)))

par(bty='n', mar=c(5,4.2,2,0))
plot_ad_seq(corr_results[[1]],672);
mtext('A',adj=0,cex=1.5)
plot_ad_seq(corr_results[[1]],672,type='late');
mtext('B',adj=0,cex=1.5)

par(mar=c(4.5,2.5,2,0))
hist(corr_nofilt$early_R, main='Assembly', ylab='', xlab=paste('Adjusted R Squared Values (n=',length(corr_nofilt$early_R_sq),')', sep=''), freq=TRUE)
mtext('C',adj=0,cex=1.5)
hist(corr_nofilt$late_R,main='Disassembly',ylab='',xlab=paste('Adjusted R Squared Values (n=',length(corr_nofilt$late_R_sq),')', sep=''), freq=TRUE)
mtext('D',adj=0,cex=1.5)

par(mar=c(2.6,5,2,0))
boxplot_with_points(list(corr_filt$early_slope,corr_filt$late_slope), names=c(paste('Assembly (n=',length(corr_filt$early_R_sq),')', sep=''), paste('Disassembly (n=',length(corr_filt$late_R_sq),')', sep='')), boxwex=0.6, ylab=expression(paste('Rate (',min^-1,')',sep='')))
mtext('E',adj=0,cex=1.5)
graphics.off()

#Spacial Data
pdf(file.path(out_folder,'spacial.pdf'))
par(bty='n',mar=c(4.2,4.1,2,0.2))
layout(rbind(c(1,2),c(3,4)))
breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,25,30,35,40);
start_props = hist(corr_filt$starting_edge_dist,xlab=expression(paste('Distance from Edge at Birth (', mu, 'm)', sep='')), main = '', ylab = '')
mtext('A',adj=0,cex=1.5)

start_props$breaks = c(start_props$breaks, max(corr_filt$ending_edge_dist))

hist(corr_filt$ending_edge_dist,xlab=expression(paste('Distance from Edge at Death (', mu, 'm)', sep='')), main = '', ylab = '')
mtext('B',adj=0,cex=1.5)

plot(corr_filt$starting_edge_dist[! is.na(corr_filt$stable_lifetime_e)],corr_filt$stable_lifetime_e[! is.na(corr_filt$stable_lifetime_e)],ylab='Stable Lifetime (min)',xlab=expression(paste('Distance from Edge at Birth (', mu, 'm)', sep='')))
mtext('C',adj=0,cex=1.5)

plot(corr_filt$ending_edge_dist[! is.na(corr_filt$stable_lifetime_l)],corr_filt$stable_lifetime_l[! is.na(corr_filt$stable_lifetime_l)],ylab='Stable Lifetime (min)',xlab=expression(paste('Distance from Edge at Death (', mu, 'm)', sep='')))
mtext('D',adj=0,cex=1.5)

graphics.off()

results_filt = filter_results(results)
corr_filt = filter_results(corr_results)
shrunk_filt = filter_results(shrunk_results)

birth_rate = list()
death_rate = list()
for (i in 1:length(results)) {
	birth_rate[[i]] = rep(0,dim(results[[i]]$exp_data)[[2]]-1);
	death_rate[[i]] = rep(0,dim(results[[i]]$exp_data)[[2]]-1);
	for (j in 1:dim(results[[i]]$exp_data)[[1]]) {
		data_line = is.na(as.vector(results[[i]]$exp_data[j,]));
		true_line = which(!data_line);
#		print(data_line)
		if (length(true_line) == 0) {
			next;
		}
		if (data_line[1]) {
			birth_rate[[i]][true_line[1]-1] = birth_rate[[i]][true_line[1]-1] + 1;
		}
		if (data_line[length(data_line)]) {
			death_rate[[i]][true_line[length(true_line)]] = death_rate[[i]][true_line[length(true_line)]] + 1;
		}		
	}
}

for (i in 1:length(birth_rate)) {
#	print(cor(birth_rate[[i]],death_rate[[i]]))
#	print(mean(birth_rate[[i]]-death_rate[[i]]))
}

for (i in 1:length(corr_filt$lin_num_early)) {
#	print(corr_filt$exp_dir[[i]]);
#	print(intersect(corr_filt$lin_num_early[[i]],corr_filt$lin_num_late[[i]]));
}