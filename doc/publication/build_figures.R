source('linear_regions.R')
library(lattice)

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
general_props = gather_general_props(results)

results_nofilt = filter_mixed_results(results, corr_results, needed_R_sq = -Inf, needed_p_val = 1.1)

results_filt = filter_mixed_results(results, corr_results)

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


########################################
#Kinetics Figure
########################################
pdf(file.path(out_folder,'kinetics.pdf'),height=9)
layout(rbind(c(1,2),c(3,4),c(5)))

par(bty='n', mar=c(5,4.2,2,0))
plot_ad_seq(corr_results[[1]],672);
mtext('A',adj=0,cex=1.5)
plot_ad_seq(corr_results[[1]],672,type='late');
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

hist(results_filt$dis$edge_dist[! is.na(results_filt$dis$stable_lifetime)] ,xlab=expression(paste('Distance from Edge at Death (', mu, 'm)', sep='')), main = '', ylab = '', breaks = start_props$breaks)
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
#Birth/Death Rates Figure
########################################

birth_rate = list()
death_rate = list()
ad_count = list()
ad_diff = list()
for (i in 1:length(results)) {
	alive_dead_mat = ! is.na(results[[i]]$exp_data)
	ad_count[[i]] = as.numeric(apply(alive_dead_mat, 2, sum))
	ad_diff[[i]] = rep(0,length(ad_count[[i]]));
	for (j in 2:length(ad_count[[i]])) {
		ad_diff[[i]][j] = ad_count[[i]][j] - ad_count[[i]][j - 1]
	}
	
	birth_rate[[i]] = rep(0,dim(results[[i]]$exp_data)[[2]]);
	death_rate[[i]] = rep(0,dim(results[[i]]$exp_data)[[2]]);
	for (j in 1:dim(results[[i]]$exp_data)[[1]]) {
		birth_time = find_birth_time(alive_dead_mat[j,])
		if (is.numeric(birth_time)) {
			birth_rate[[i]][birth_time] = birth_rate[[i]][birth_time] + 1
		}
		death_time = find_death_time(alive_dead_mat[j,])
		if (is.numeric(death_time)) {
			death_rate[[i]][death_time] = death_rate[[i]][death_time] + 1
		}		
	}
}

corrs = list()
corrs$x = 1:length(birth_rate)
for (i in 1:length(birth_rate)) {
	temp = cor.test(birth_rate[[i]], death_rate[[i]])
	corrs$y = c(corrs$y, as.numeric(temp$est))
	corrs$yplus = c(corrs$yplus, temp$conf.int[2])
	corrs$yminus = c(corrs$yminus, temp$conf.int[1])
}

library(Hmisc)
pdf(file.path(out_folder,'birth_death_corr.pdf'))
par(bty='n',mar=c(4.2,4.1,2,0.2))
errbar(corrs$x,corrs$y,corrs$yplus,corrs$yminus, ylab ='Correlation', xlab ='Exp Number')
graphics.off()