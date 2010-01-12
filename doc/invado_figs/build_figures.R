rm(list = ls())
source('../../src_invado/analyze_cell_features/invado_analysis_lib.R')
source('../../src_invado/analyze_cell_features/FA_analysis_lib.R')
library(lattice)
library(geneplotter)
library(Hmisc)
library(boot)

debug=TRUE;

################################################################################
# Data Reading/Initial Processing
################################################################################

control_props = gather_invado_properties(Sys.glob('../../results/Invadopodia/control/*/adhesion_props/'));
coro_kd_props = gather_invado_properties(Sys.glob('../../results/Invadopodia/Coro1C/*/adhesion_props/'));

subsets = list();

subsets$invado_control = subset(control_props, high_conf_int < 0);
subsets$invado_coro_kd = subset(coro_kd_props, high_conf_int < 0);
subsets$non_invado_control = subset(control_props, high_conf_int >= 0);
subsets$non_invado_coro_kd = subset(coro_kd_props, high_conf_int >= 0);

stop()

################################################################################
# Plotting
################################################################################

########################################
# Longevity Bi-hist
########################################
bin_size = 50;
max_data = 5*max(subsets$invado_control$longevity, subsets$invado_coro_kd$longevity);
max_data = ceil(max_data/bin_size)*bin_size;
min_data = 5*min(subsets$invado_control$longevity, subsets$invado_coro_kd$longevity);

data_1 = hist(subsets$invado_control$longevity*5, breaks=seq(min_data,max_data,by=bin_size));
data_2 = hist(subsets$invado_coro_kd$longevity*5, breaks=data_1$breaks);
graphics.off()

all_data = rbind(data_1$counts/sum(data_1$counts), data_2$counts/sum(data_2$counts));
range_names = c()
for (i in 1:(length(data_1$breaks)-1)) {
    range_names = c(range_names, paste(data_1$breaks[i],"-",data_1$breaks[i+1],sep=''));
}

svg('longevity_hist.svg')
par(mar = c(5.75,4,0.25,0))
barplot(all_data, beside=T, names = range_names, ylab="Percentage", 
        #las controls the orientation of the labels, las=2 is perpendicular to axis
        legend = c('Control', 'Coronin1C KD'), las=2)
mtext("Longevity Ranges (min)", side=1, line=4.5)
graphics.off()

plot(subsets$invado_control$longevity*5, subsets$invado_control$largest_area, pch=20, cex=0.75)
points(subsets$invado_coro_kd$longevity*5, subsets$invado_coro_kd$largest_area, col='red', pch=20, cex=0.75)

t.test(subsets$invado_control$longevity*5, subsets$invado_coro_kd$longevity*5)
t.test(subsets$invado_control$largest_area, subsets$invado_coro_kd$largest_area)

########################################
# Largest Size Bi-hist
########################################
bin_size = 0.5;
max_data = max(subsets$invado_control$largest_area, subsets$invado_coro_kd$largest_area);
max_data = ceil(max_data/bin_size)*bin_size;
min_data = min(subsets$invado_control$largest_area, subsets$invado_coro_kd$largest_area);
min_data = floor(min_data/bin_size)*bin_size;

data_1 = hist(subsets$invado_control$largest_area, breaks=seq(min_data,max_data,by=bin_size));
data_2 = hist(subsets$invado_coro_kd$largest_area, breaks=data_1$breaks);
graphics.off()

all_data = rbind(data_1$counts/sum(data_1$counts), data_2$counts/sum(data_2$counts));
range_names = c()
for (i in 1:(length(data_1$breaks)-1)) {
    range_names = c(range_names, paste(data_1$breaks[i],"-",data_1$breaks[i+1],sep=''));
}

svg('largest_area_hist.svg')
par(mar = c(4.75,4,0.25,0))
barplot(all_data, beside=T, names = range_names, ylab="Percentage", 
        #las controls the orientation of the labels, las=2 is perpendicular to axis
        legend = c('Control', 'Coronin1C KD'), las=2)
mtext('Max Invadopodia Area (\u03BCm\u00B2)', side=1, line=3.5)
graphics.off()

########################################
# Average Local Diff Values Bi-hist
########################################
bin_size = 0.5E-3*10^2;

control_mean_vals = subsets$invado_control$mean_vals*-10^2;
coro_kd_mean_vals = subsets$invado_coro_kd$mean_vals*-10^2;

max_data = max(control_mean_vals, coro_kd_mean_vals);
max_data = ceil(max_data/bin_size)*bin_size;
min_data = min(control_mean_vals, coro_kd_mean_vals);
min_data = floor(min_data/bin_size)*bin_size;

data_1 = hist(control_mean_vals, breaks=seq(min_data,max_data,by=bin_size));
data_2 = hist(coro_kd_mean_vals, breaks=data_1$breaks);
graphics.off()

all_data = rbind(data_1$counts/sum(data_1$counts), data_2$counts/sum(data_2$counts));
range_names = c()
for (i in 1:(length(data_1$breaks)-1)) {
    range_names = c(range_names, paste(data_1$breaks[i],"-",data_1$breaks[i+1],sep=''));
}

svg('mean_vals_hist.svg')
par(mar = c(6,4,0.25,0))
barplot(all_data, beside=T, names = range_names, ylab="Percentage", 
        #las controls the orientation of the labels, las=2 is perpendicular to axis
        las=2)
mtext('Mean Local Difference (arbitrary units * -10\u00B2)', side=1, line=4.75)
legend('topright', c('Control', 'Coronin1C KD'), fill = c('gray30','gray90'))
graphics.off()

plot(subsets$invado_control$longevity, subsets$invado_control$mean_vals, ylim=c(min_data,max_data))
points(subsets$invado_coro_kd$longevity, subsets$invado_coro_kd$mean_vals, col='red')

