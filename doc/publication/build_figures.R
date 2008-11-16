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
hist(corr_nofilt$early$R, main='Assembly', ylab='', xlab=paste('Adjusted R Squared Values (n=',length(corr_nofilt$early_R_sq),')', sep=''), freq=TRUE)
mtext('C',adj=0,cex=1.5)
hist(corr_nofilt$late$R,main='Disassembly',ylab='',xlab=paste('Adjusted R Squared Values (n=',length(corr_nofilt$late_R_sq),')', sep=''), freq=TRUE)
mtext('D',adj=0,cex=1.5)

par(mar=c(2.6,5,2,0))
boxplot_with_points(list(corr_filt$early$slope,corr_filt$late$slope), names=c(paste('Assembly (n=',length(corr_filt$early$R_sq),')', sep=''), paste('Disassembly (n=',length(corr_filt$late$R_sq),')', sep='')), boxwex=0.6, ylab=expression(paste('Rate (',min^-1,')',sep='')))
mtext('E',adj=0,cex=1.5)
graphics.off()

#Spacial Data
pdf(file.path(out_folder,'spacial.pdf'))
par(bty='n',mar=c(4.2,4.1,2,0.2))
layout(rbind(c(1,2),c(3,4)))
breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,25,30,35,40);
start_props = hist(corr_filt$early$edge_dist,xlab=expression(paste('Distance from Edge at Birth (', mu, 'm)', sep='')), main = '', ylab = '')
mtext('A',adj=0,cex=1.5)

start_props$breaks = c(start_props$breaks, max(corr_filt$late$edge_dist))

hist(corr_filt$late$edge_dist,xlab=expression(paste('Distance from Edge at Death (', mu, 'm)', sep='')), main = '', ylab = '')
mtext('B',adj=0,cex=1.5)

plot(corr_filt$early$edge_dist[! is.na(corr_filt$early$stable_lifetime)],corr_filt$early$stable_lifetime[! is.na(corr_filt$early$stable_lifetime)],ylab='Stable Lifetime (min)',xlab=expression(paste('Distance from Edge at Birth (', mu, 'm)', sep='')))
mtext('C',adj=0,cex=1.5)

plot(corr_filt$late$edge_dist[! is.na(corr_filt$late$stable_lifetime)],corr_filt$late$stable_lifetime[! is.na(corr_filt$late$stable_lifetime)],ylab='Stable Lifetime (min)',xlab=expression(paste('Distance from Edge at Death (', mu, 'm)', sep='')))
mtext('D',adj=0,cex=1.5)

graphics.off()