rm(list = ls())
source('FA_analysis_lib.R')
library(lattice)
library(geneplotter)
library(Hmisc)
debug=TRUE;

################################################################################
#Result loading
################################################################################
raw_data <- list()
single_props <- list()

#Wild-type FA
exp_dirs <- Sys.glob('../../results/focal_adhesions/*/adhesion_props/models/')
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]

raw_data$wild_type$intensity = load_results(exp_dirs,file.path('intensity.Rdata'));
raw_data$wild_type$corrected_intensity = load_results(exp_dirs,file.path('local_corrected.Rdata'));
raw_data$wild_type$static_props <- load_data_files(exp_dirs, 
    file.path('..','individual_adhesions.csv'), headers=T, debug=FALSE, inc_exp_names=FALSE);

# #S178A Results
exp_dirs_S <- Sys.glob('../../results/S178A/*/adhesion_props/models/')
exp_dirs_S <- exp_dirs_S[file_test('-d',exp_dirs_S)]

raw_data$S178A$intensity = load_results(exp_dirs_S,file.path('intensity.Rdata'));
raw_data$S178A$corrected_intensity = load_results(exp_dirs_S,file.path('local_corrected.Rdata'));
raw_data$S178A$static_props <- load_data_files(exp_dirs_S, 
    file.path('..','individual_adhesions.csv'), headers=T, debug=FALSE, inc_exp_names=FALSE);

print('Done Loading Data')

########################################
#Result filtering
########################################

processed = list();
dynamic_props = list();
static_props = list();
for (exp_type in names(raw_data)) {
    for (property in names(raw_data[[exp_type]])) {
        if (property == "static_props") {
            next;
        }
        if (debug) {
            print(paste("Filtering", exp_type, property));
        }

        processed$no_filt[[exp_type]][[property]] = filter_results(raw_data[[exp_type]][[property]], 
            min_R_sq = -Inf, max_p_val = Inf);
        
        processed$only_signif[[exp_type]][[property]] = filter_results(raw_data[[exp_type]][[property]], 
            min_R_sq = -Inf, max_p_val = 0.05);
        
        processed$high_Rsq[[exp_type]][[property]] = filter_results(raw_data[[exp_type]][[property]]);
    }
    
    dynamic_props[[exp_type]] = gather_general_dynamic_props(raw_data[[exp_type]]$intensity)
    static_props[[exp_type]] = gather_static_props(raw_data[[exp_type]]$static_props) 
}

print('Done Filtering Data')
out_folder = '../../doc/publication/figures'
dir.create(out_folder,recursive=TRUE, showWarnings=FALSE);

stop()

################################################################################
#Plotting
################################################################################
out_folder = '../../doc/publication/figures'
dir.create(out_folder,recursive=TRUE, showWarnings=FALSE);


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
area_data = static_props$wild_type$Area[static_props$wild_type$Area < 5];
area_hist = hist(area_data, main="", ylab = "FA Count", 
	 xlab = expression(paste('FA Area (', symbol("m"), m^2, ')',sep='')));
    
#area_hist_model = lm(log(area_hist$counts) ~area_hist$mids);
#predictions = exp(predict(area_hist_model))
#lines(area_hist$mids,predictions, col='red', lwd=2)

mtext('D',adj=-.2,side=3,line=-0.25,cex=1.5)
hist(static_props$wild_type$ad_sig, main="", ylab = "FA Count", xlab = "Normalized Average Paxillin Intensity")

par(bty='n', mar=c(4,4.2,1.2,0.1))
hist(static_props$wild_type$ax[static_props$wild_type$ax < 8], main="", ylab = "FA Count",  xlab = "Axial Ratio")
hist(static_props$wild_type$cent_dist, main="", ylab = "FA Count",  
	 xlab = expression(paste("Distance from Edge (", symbol("m"), m, ')',sep='')))
hist(dynamic_props$wild_type$longevity, main="", ylab = "FA Count",  xlab = "Longevity (min)")
graphics.off()

svg(file.path(out_folder,'statics','longev_inset.svg'), width=3, height=3/2)
par(bty='n', mar=c(2,2,0.5,0))
hist(dynamic_props$wild_type$longevity[!is.na(dynamic_props$wild_type$longevity) & 
                                        dynamic_props$wild_type$longevity > 20], 
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

raw_data_wt_one = raw_data$wild_type$intensity[[1]]

par(bty='n', mar=c(4,4.2,1.5,0))

plot.new()
mtext('A',adj=-.19,side=3,line=-0.5,cex=1.5)

ad_num = 1799
ad_num = 675
plot_ad_seq(raw_data_wt_one, ad_num, type='overall', 
	phase_lengths=c(raw_data_wt_one$assembly$length[ad_num],raw_data_wt_one$disassembly$length[ad_num]))
mtext('B',adj=-.19,side=3,line=-0.5,cex=1.5)

par(bty='n', mar=c(4,4.2,4,0))
plot_ad_seq(raw_data_wt_one, ad_num, main = 'Assembly');
limits = par("usr");
text(3,(limits[4]-limits[3])*0.8+limits[3],pos=3,expression(paste(R^2,' = 0.920')))
text(3,(limits[4]-limits[3])*0.8+limits[3],pos=3, offset=c(-0.65,0),
	substitute(paste('Slope = ', x), list(x=sprintf('%.03f',raw_data_wt_one$assembly$slope[ad_num]))))
mtext('C',adj=-.19,side=3,line=-0.5,cex=1.5)

plot_ad_seq(raw_data_wt_one,ad_num,type='disassembly', main = 'Disassembly')
limits = par("usr");
text(3,(limits[4]-limits[3])*0.8+limits[3],pos=3, expression(paste(R^2,' = 0.961')))
text(3,(limits[4]-limits[3])*0.8+limits[3], offset=c(-0.65,0),
	substitute(paste('Slope = ', x), list(x=sprintf('%.03f',raw_data_wt_one$disassembly$slope[ad_num]))))
mtext('D',adj=-.19,side=3,line=-0.5,cex=1.5)
par(bty='n', mar=c(2.1,4.2,1.1,0))

boxplot_with_points(list(processed$high_Rsq$wild_type$intensity$assembly$slope,
                         processed$high_Rsq$wild_type$intensity$disassembly$slope),
		    names=c('Assembly', 'Disassembly'), boxwex=0.6, 
		    ylab=expression(paste('Rate (',min^-1,')',sep='')), point_cex=0.6, with.median.props=FALSE)
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

R_sq_data = processed$only_signif$wild_type$intensity;

par(bty='n', mar=c(4,4.2,2,0))

hist(R_sq_data$assembly$R_sq, main='Assembly', freq=TRUE,
	 xlab=paste('Adjusted R Squared Values (n=',length(R_sq_data$assembly$R_sq),')', sep=''), 
	 ylab='# of Focal Adhesions')

plot_dims = par("usr");
sorted_r_vals = sort(R_sq_data$assem$R)
segments(sorted_r_vals[floor(length(sorted_r_vals)/2)],0, 
         sorted_r_vals[floor(length(sorted_r_vals)/2)],plot_dims[4], col='red', lwd=2)
mtext('A',adj=-.2,side=3,line=-0.5,cex=1.5)

hist(R_sq_data$dis$R,main='Disassembly', freq=TRUE,
	 xlab=paste('Adjusted R Squared Values (n=',length(R_sq_data$d$R_sq),')', sep=''), 
	 ylab='# of Focal Adhesions')

plot_dims = par("usr");
sorted_r_vals = sort(R_sq_data$dis$R)
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

wt_high_Rsq = processed$high_Rsq$wild_type$intensity;

breaks_end = ceil(max(c(wt_high_Rsq$a$edge_dist,wt_high_Rsq$dis$edge_dist), na.rm=T));
if (breaks_end %% 2 != 0) {
	breaks_end = breaks_end + 1;
}
these_breaks = seq(0,breaks_end,by=2);

hist(wt_high_Rsq$a$edge_dist, 
     xlab=expression(paste('Distance from Edge at Birth (', mu, 'm) n=309', sep='')), 
     main = '', ylab = '# of Focal Adhesions', breaks=these_breaks)
mtext('A',adj=-.25,side=3,line=-1.5,cex=1.5)

par(bty='n',mar=c(4.2,4.1,0.1,0))
hist(wt_high_Rsq$dis$edge_dist,
     xlab=expression(paste('Distance from Edge at Death (', mu, 'm) n=470', sep='')), 
     main = '', ylab = '# of Focal Adhesions', breaks=these_breaks)
mtext('B',adj=-.25,side=3,line=-1.5,cex=1.5)

par(bty='n',mar=c(4.2,4.1,0.1,1))
max_rate = max(c(wt_high_Rsq$a$slope, wt_high_Rsq$d$slope));
plot(wt_high_Rsq$a$edge_dist, pch=19, cex=0.5,
	 wt_high_Rsq$a$slope,
	 xlim = c(0,breaks_end),
         ylim = c(0,max_rate),
	 ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')),
	 xlab=expression(paste('Distance from Edge at Birth (', mu, 'm)', sep='')))
mtext('C',adj=-.25,side=3,line=-1.5,cex=1.5)

par(bty='n',mar=c(4.2,4.1,0.1,0))
plot(wt_high_Rsq$d$edge_dist, pch=19, cex = 0.5,
	 wt_high_Rsq$d$slope, 
	 xlim = c(0,breaks_end),
         ylim = c(0,max_rate),
	 ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')),
	 xlab=expression(paste('Distance from Edge at Death (', mu, 'm)', sep='')))
mtext('D',adj=-.25,side=3,line=-1.5,cex=1.5)
graphics.off()

#Alternate combined plots
dir.create(dirname(file.path(out_folder,'spatial','spatial_alt.svg')), 
    recursive=TRUE, showWarnings=FALSE);
svg(file.path(out_folder,'spatial','spatial_alt.svg'),height=4, width=10)
par(bty='n',mar=c(4.2,4.1,0.1,6))
layout(rbind(c(1,2)))

breaks_end = ceil(max(c(wt_high_Rsq$a$edge_dist,wt_high_Rsq$dis$edge_dist), na.rm=T));
if (breaks_end %% 2 != 0) {
	breaks_end = breaks_end + 1;
}
these_breaks = seq(0,breaks_end,by=2);
max_rate = max(c(wt_high_Rsq$a$slope, wt_high_Rsq$d$slope));

#assembly rates
hist(wt_high_Rsq$a$edge_dist, 
     xlab=expression(paste('Distance from Edge at Birth (', mu, 'm) n=309', sep='')), 
     main = '', ylab = '# of Focal Adhesions', breaks=these_breaks)
plot_range = par('usr')
# segments(1.5,0,1.5,plot_range[4], col='purple', lwd=2)
points(wt_high_Rsq$a$edge_dist, pch=19, cex=0.35, col='darkgreen',
	 wt_high_Rsq$a$slope*((plot_range[4]-10)/max_rate),
	 ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')))
axis_ticks = axTicks(2);
axis(4, at = axis_ticks*((plot_range[4]-10)/max_rate), labels=axis_ticks, col='darkgreen')
mtext(expression(paste('Assembly Rate (',min^-1,')',sep='')),side=4,line=3, col='darkgreen');
mtext('A',adj=-.22,side=3,line=-1.5,cex=1.75)

#disassembly rates
par(bty='n',mar=c(4.2,4.1,0.1,4))
hist(wt_high_Rsq$dis$edge_dist,
     xlab=expression(paste('Distance from Edge at Death (', mu, 'm) n=470', sep='')), 
     main = '', ylab = '# of Focal Adhesions', breaks=these_breaks)
plot_range = par('usr')
# segments(1.5,0,1.5,plot_range[4], col='purple', lwd=2)
points(wt_high_Rsq$d$edge_dist, pch=19, cex=0.35, col='red',
	 wt_high_Rsq$d$slope*((plot_range[4]-10)/max_rate),
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

wt_no_filt = processed$no_filt$wild_type$intensity;

breaks_end = ceil(max(c(wt_no_filt$a$edge_dist,wt_no_filt$dis$edge_dist), na.rm=T));
if (breaks_end %% 2 != 0) {
	breaks_end = breaks_end + 1;
}
these_breaks = seq(0,breaks_end,by=2);

hist(wt_no_filt$a$edge_dist, 
    xlab=expression(paste('Distance from Edge at Birth (', mu, 'm) n=1169', sep='')),
    main = '', ylab = '# of Focal Adhesions', breaks=these_breaks)
mtext('A',adj=-.2,side=3,line=-0.5,cex=1.5)

hist(wt_no_filt$dis$edge_dist,
     xlab=expression(paste('Distance from Edge at Death (', mu, 'm) n=1460', sep='')), 
     main = '', ylab = '# of Focal Adhesions',breaks=these_breaks)
mtext('B',adj=-.2,side=3,line=-0.5,cex=1.5)

par(bty='n',mar=c(4.2,4.1,2,0.4))

smoothScatter(wt_no_filt$a$edge_dist, 
	 wt_no_filt$a$slope,
	 xlim = c(0,breaks_end), 
	 ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')),
	 xlab=expression(paste('Distance from Edge at Birth (', mu, 'm)', sep='')))
mtext('C',adj=-.2,side=3,line=-0.5,cex=1.5)

par(bty='n',mar=c(4.2,4.1,2,0.95))

smoothScatter(wt_no_filt$d$edge_dist, 
	 wt_no_filt$d$slope,
	 xlim = c(0,breaks_end), 
	 ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')),
	 xlab=expression(paste('Distance from Edge at Death (', mu, 'm)', sep='')))
mtext('D',adj=-.2,side=3,line=-0.5,cex=1.5)
graphics.off()


svg(file.path(out_folder,'supplemental','birth_vs_death_pos.svg'))
par(bty='n',mar=c(4.2,4.1,2,0.2))


plot(wt_no_filt$j$birth_dist, wt_no_filt$j$death_dist, 
     xlab=expression(paste('Distance from Edge at Birth (', mu, 'm)', sep='')),
     ylab=expression(paste('Distance from Edge at Death (', mu, 'm)', sep='')), 
     pch=20, cex=0.75)
#points(results_S_nofilt$j$b, results_S_nofilt$j$d, col='red', pch=20, cex=0.5)

birth_vs_death_model <- lm(death_dist ~ birth_dist, data=wt_no_filt$joint)
abline(birth_vs_death_model, col='green', lwd = 3)
model_summary <- summary(birth_vs_death_model)
x_data <- data.frame(birth_dist = seq(min(wt_no_filt$j$b, na.rm=T),max(wt_no_filt$j$b, na.rm=T),by=0.01))
line_conf = predict(birth_vs_death_model, x_data, interval="confidence", level=0.95)
lines(x_data$birth_dist, line_conf[,2], col='red', lty=2, lwd = 3)
lines(x_data$birth_dist, line_conf[,3], col='red', lty=2, lwd = 3)

segments(0,0,max(wt_no_filt$j$b, na.rm=T), max(wt_no_filt$j$b, na.rm=T), col='blue', lty=4, lwd = 3)

graphics.off()

svg(file.path(out_folder,'supplemental','birth_vs_death_pos_S178A.svg'))
par(bty='n',mar=c(4.2,4.1,2,0.2))

S178A_no_filt = processed$no_filt$S178A$intensity;

plot(S178A_no_filt$j$birth_dist, S178A_no_filt$j$death_dist, 
     xlab=expression(paste('Distance from Edge at Birth (', mu, 'm)', sep='')),
     ylab=expression(paste('Distance from Edge at Death (', mu, 'm)', sep='')), 
     pch=20, cex=0.75)
#points(S178A_no_filt$j$b, S178A_no_filt$j$d, col='red', pch=20, cex=0.5)

birth_vs_death_model <- lm(death_dist ~ birth_dist, data=S178A_no_filt$joint)
abline(birth_vs_death_model, col='green', lwd = 3)
model_summary <- summary(birth_vs_death_model)
x_data <- data.frame(birth_dist = seq(min(S178A_no_filt$j$b, na.rm=T),max(S178A_no_filt$j$b, na.rm=T),by=0.01))
line_conf = predict(birth_vs_death_model, x_data, interval="confidence", level=0.95)
lines(x_data$birth_dist, line_conf[,2], col='red', lty=2, lwd = 3)
lines(x_data$birth_dist, line_conf[,3], col='red', lty=2, lwd = 3)

segments(0,0,max(S178A_no_filt$j$b, na.rm=T), max(S178A_no_filt$j$b, na.rm=T), col='blue', lty=4, lwd = 3)

graphics.off()
print('Done with Spacial')

############################################################
#Comparing S178A to Wild-type
############################################################
dir.create(dirname(file.path(out_folder,'lifetimes','stable_mean.pdf')), 
    recursive=TRUE, showWarnings=FALSE);
pdf(file.path(out_folder,'lifetimes','stable_mean.pdf'), width=8.5, pointsize=8)
boxplot_with_points(
        list(processed$only_signif$wild_type$corrected_intensity$joint$stable_mean, 
             processed$only_signif$S178A$corrected_intensity$joint$stable_mean, 
             processed$only_signif$wild_type$intensity$joint$stable_mean, 
             processed$only_signif$S178A$intensity$joint$stable_mean
            ),
        names=c('Local WT','Local S189A', 'WT','S178A'), ylab='Mean Stable Intensity', with.median.props=FALSE)
graphics.off()

print('Done with Stable Lifetime Averages')

########################################
#Lifetime Phases
########################################

stage_data <- gather_stage_lengths(processed$only_signif$wild_type$corrected_intensity, 
                                   processed$only_signif$S178A$corrected_intensity)

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
plot_signif_bracket(upper_left,lower_right, orientation='upside_down', over_text = "**", text_x_adj=-0.005)

#right side up identification bar
upper_left = c(err_bars[3,1], max(err_bars[3,4],err_bars[6,4]) + sep_from_data + bar_length);
lower_right = c(err_bars[6,1], max(err_bars[3,4],err_bars[6,4]) + sep_from_data);
plot_signif_bracket(upper_left,lower_right, over_text = "*", text_x_adj=-0.005)

graphics.off()

svg(file.path(out_folder,'lifetimes','adhesion_phase_lifetimes_alt.svg'))
par(bty='n',mar=c(2.1,4,0,0))

sideways_err_bars = plot_stage_length_data(stage_data, type='side_by_side', top_gap=1.5, 
    names=c('Wild-type','S178A'), sideways_high_xlim=11);

#Signifcance Bars
bar_length = 0.4;
sep_from_data = 0.2;

#assembly significance bar
upper_left = c(sideways_err_bars[1,1], max(sideways_err_bars[1,4],sideways_err_bars[4,4]) + sep_from_data + bar_length);
lower_right = c(sideways_err_bars[4,1], max(sideways_err_bars[1,4],sideways_err_bars[4,4]) + sep_from_data);
plot_signif_bracket(upper_left,lower_right, over_text = "**")

#disassembly significance bar
upper_left = c(sideways_err_bars[3,1], max(sideways_err_bars[3,4],sideways_err_bars[6,4]) + sep_from_data + bar_length);
lower_right = c(sideways_err_bars[6,1], max(sideways_err_bars[3,4],sideways_err_bars[6,4]) + sep_from_data);
plot_signif_bracket(upper_left,lower_right, over_text = "*")

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

wt_high_Rsq = processed$high_Rsq$wild_type$intensity;
S178A_high_Rsq = processed$high_Rsq$S178A$intensity;

max_rate = max(c(wt_high_Rsq$as$slope,S178A_high_Rsq$as$slope,
                 wt_high_Rsq$dis$slope,S178A_high_Rsq$dis$slope));

#Panel Assembly Rates
boxplot_with_points(list(wt_high_Rsq$as$slope,S178A_high_Rsq$as$slope), 
        names=c('Wild-type','S178A'), 
        colors=c('orange','blue'),
        ylim = c(0,max_rate + 0.012),
        ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')),
        median.props.color='brown'
)

# bar_length = .005;
# sep_from_data = 0.005;
# 
# upper_left = c(1, max_rate + sep_from_data + bar_length);
# lower_right = c(2, max_rate + sep_from_data);
# lines(c(upper_left[1],upper_left[1],lower_right[1], lower_right[1]),
# 	  c(lower_right[2],upper_left[2],upper_left[2],lower_right[2]))  
# text(mean(c(upper_left[1],lower_right[1]))-0.005,upper_left[2]+sep_from_data,"***",cex=1.5)
mtext('A',adj=-.25,side=3,line=-1.5,cex=1.5)	    

#Panel Disassembly Rates
boxplot_with_points(list(wt_high_Rsq$dis$slope,S178A_high_Rsq$dis$slope), 
    names=c('Wild-type','S178A'), 
    colors=c('orange','blue'),
    ylim = c(0,max_rate + 0.012),
    ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')),
    median.props.color='brown'

)

# upper_left = c(1, max_rate + sep_from_data + bar_length);
# lower_right = c(2, max_rate + sep_from_data);
# lines(c(upper_left[1],upper_left[1],lower_right[1], lower_right[1]),
# 	  c(lower_right[2],upper_left[2],upper_left[2],lower_right[2]))  
# text(mean(c(upper_left[1],lower_right[1]))-0.005,upper_left[2]+sep_from_data,"*",cex=1.5)
mtext('B',adj=-.25,side=3,line=-1.5,cex=1.5)

#Panel Birth Distances
par(bty='n', mar=c(2,4,1.5,0))
max_dist = max(c(wt_high_Rsq$as$edge_dist,S178A_high_Rsq$as$edge_dist,
                 wt_high_Rsq$dis$edge_dist,S178A_high_Rsq$dis$edge_dist), na.rm=T)
boxplot_with_points(list(wt_high_Rsq$as$edge_dist,S178A_high_Rsq$as$edge_dist), 
    names=c('Wild-type','S178A'), 
    ylim=c(0,max_dist+2), 
    colors=c('orange','blue'),
    ylab=expression(paste('Distance from Edge at Birth (',mu,'m)',sep='')),
    median.props.color='brown'
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
boxplot_with_points(list(wt_high_Rsq$dis$edge_dist,S178A_high_Rsq$dis$edge_dist), 
    names=c('Wild-type','S178A'), 
    ylim=c(0,max_dist+2),
    colors=c('orange','blue'),
    ylab=expression(paste('Distance from Edge at Death (',mu,'m)',sep='')),
    median.props.color='brown'
)
mtext('D',adj=-.25,side=3,line=-1.5,cex=1.5)	    

graphics.off()

###############
#Supplemental
###############
dir.create(dirname(file.path(out_folder,'supplemental','unfilt_S178A_vs_wild-type.pdf')), 
    recursive=TRUE, showWarnings=FALSE);
pdf(file.path(out_folder,'supplemental','unfilt_S178A_vs_wild-type.pdf'))
layout(rbind(c(1,2),c(3,4)))
par(bty='n', mar=c(2,4,0,0))

wt_only_signif = processed$only_signif$wild_type$intensity;
S178A_only_signif = processed$only_signif$S178A$intensity;

max_rate = max(c(wt_only_signif$as$slope,S178A_only_signif$as$slope,
                 wt_only_signif$dis$slope,S178A_only_signif$dis$slope));

#Panel Assembly Rates
boxplot_with_points(list(wt_only_signif$as$slope,S178A_only_signif$as$slope), 
    names=c('Wild-type','S178A'), 
    colors=c('orange','blue'),
    ylim = c(0,max_rate),
    ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')),
    inc.points=FALSE,
    median.props.pos = rbind(c(0.5,0.5),c(0,0.9))
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
boxplot_with_points(list(wt_only_signif$dis$slope,S178A_only_signif$dis$slope), 
    names=c('Wild-type','S178A'), 
    colors=c('orange','blue'),
    ylim = c(0,max_rate),
    ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')),
    inc.points = FALSE,
    median.props.pos = rbind(c(0.5,0.5),c(0,0.9))
)

#upper_left = c(1, max_rate + sep_from_data + bar_length);
#lower_right = c(2, max_rate + sep_from_data);
#lines(c(upper_left[1],upper_left[1],lower_right[1], lower_right[1]),
#	  c(lower_right[2],upper_left[2],upper_left[2],lower_right[2]))  
#text(mean(c(upper_left[1],lower_right[1]))-0.005,upper_left[2]+sep_from_data,"*",cex=1.5)
mtext('B',adj=-.25,side=3,line=-1.5,cex=1.5)

#Panel Birth Distances
par(bty='n', mar=c(2,4,1.5,0))
max_dist = max(c(wt_only_signif$as$edge_dist,S178A_only_signif$as$edge_dist,
                 wt_only_signif$dis$edge_dist,S178A_only_signif$dis$edge_dist), na.rm=T)
boxplot_with_points(list(wt_only_signif$as$edge_dist,S178A_only_signif$as$edge_dist), 
    names=c('Wild-type','S178A'), 
    ylim=c(0,max_dist), 
    colors=c('orange','blue'),
    ylab=expression(paste('Distance from Edge at Birth (',mu,'m)',sep='')),
    inc.points = FALSE,
    median.props.pos = rbind(c(0.5,0.5),c(0,0.9))
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
boxplot_with_points(list(wt_only_signif$dis$edge_dist,S178A_only_signif$dis$edge_dist), 
    names=c('Wild-type','S178A'), 
    ylim=c(0,max_dist),
    colors=c('orange','blue'),
    ylab=expression(paste('Distance from Edge at Death (',mu,'m)',sep='')), 
    inc.points = FALSE,
    median.props.pos = rbind(c(0.5,0.5),c(0,0.9))
)
mtext('D',adj=-.25,side=3,line=-1.5,cex=1.5)	    

graphics.off()
print('Done with S178A Comparisons')

############################################################
#Distance versus Background corrected correlation
############################################################

# all = list()
# for (i in 1:length(background_correlation)) {
# 	all$distances = c(all$distances, background_correlation[[i]]$distances)
# 	all$correlations = c(all$correlations, background_correlation[[i]]$correlations)
# }	
# all = as.data.frame(all)
# 
# all_S = list()
# for (i in 1:length(background_correlation_S)) {
# 	all_S$distances = c(all_S$distances, background_correlation_S[[i]]$distances)
# 	all_S$correlations = c(all_S$correlations, background_correlation_S[[i]]$correlations)
# }	
# all_S = as.data.frame(all_S)
# 
# binned_all = bin_corr_data(all, bin_max=200)
# binned_all_S = bin_corr_data(all_S, bin_max=200)
# 
# dir.create(dirname(file.path(out_folder,'S178A','dist_vs_corr.pdf')), 
#     recursive=TRUE, showWarnings=FALSE);
# pdf(file.path(out_folder,'S178A','dist_vs_corr.pdf'))
# y_max = max(c(binned_all$upper, binned_all_S$upper), na.rm=TRUE)
# 
# plot(binned_all$mids, binned_all$means, ylim=c(0, y_max), 
#         xlab=expression(paste('Mean Distance (', mu, 'm)', sep='')), ylab="correlation")
# errbar(binned_all$mids, binned_all$means, binned_all$upper, binned_all$lower, add=TRUE, xlab="", ylab="")
# 
# #plot(binned_all_S$mids, binned_all_S$means, ylim=c(0, y_max), col='red', add=TRUE)
# errbar(binned_all_S$mids, binned_all_S$means, binned_all_S$upper, binned_all_S$lower, add=TRUE,col='red', xlab="", ylab="")
# graphics.off()
# 
# print('Done with Distance versus Pax Concentration Correlation')
