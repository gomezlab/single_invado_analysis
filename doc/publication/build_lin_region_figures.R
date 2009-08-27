rm(list = ls())
source('FA_analysis_lib.R')
library(lattice)
library(geneplotter)
library(Hmisc)

################################################################################
#Result loading
################################################################################
raw_data <- list();

for (i in 5:9) {
	exp_dirs <- Sys.glob(file.path('../../results/lin_region*/FA/',i,'*/adhesion_props/models/'));
	exp_dirs <- exp_dirs[file_test('-d',exp_dirs)];
	raw_data$wild[[paste("results_",i,sep='')]]  = load_results(exp_dirs,file.path('intensity.Rdata'));
	
	exp_dirs <- Sys.glob(file.path('../../results/lin_region*/S178A/',i,'*/adhesion_props/models/'));
	exp_dirs <- exp_dirs[file_test('-d',exp_dirs)];
	raw_data$S178A[[paste("results_",i,sep='')]]  = load_results(exp_dirs,file.path('intensity.Rdata'));
}

exp_dirs <- Sys.glob(file.path('../../results/focal_adhesions/*/adhesion_props/models/'));
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)];
raw_data$wild[[paste("results_",10,sep='')]]  = load_results(exp_dirs,file.path('intensity.Rdata'));

exp_dirs <- Sys.glob(file.path('../../results/S178A/*/adhesion_props/models/'));
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)];
raw_data$S178A[[paste("results_",10,sep='')]]  = load_results(exp_dirs,file.path('intensity.Rdata'));

print('Done Loading Data')

########################################
#Result filtering
########################################
processed = list();

for (exp_type in names(raw_data)) {

    for (n in names(raw_data$wild)) {
    	processed$no_filt[[exp_type]][[n]] = filter_results(raw_data[[exp_type]][[n]], 
			min_R_sq = -Inf, max_p_val = Inf)
	    
	    #now figure out the minimum length of the experiment and check to make
	    #sure that is indeed the min
	    regex_hit = regexpr('[[:digit:]]+$',n);
	    min_length = substr(n, regex_hit[1], regex_hit[1] + attr(regex_hit,"match.length"));

	    stopifnot(min(processed$no_filt[[exp_type]][[n]]$assembly$length) == min_length);
	    stopifnot(min(processed$no_filt[[exp_type]][[n]]$disassembly$length) == min_length);
	
		processed$only_signif[[exp_type]][[n]] = filter_results(raw_data[[exp_type]][[n]], 
			min_R_sq = -Inf, max_p_val = 0.05);
		
		processed$high_Rsq[[exp_type]][[n]] = filter_results(raw_data[[exp_type]][[n]]);
    }
}

print('Done Filtering Data')

################################################################################
#Plotting
################################################################################
out_folder = '../../doc/publication/figures'
dir.create(out_folder,recursive=TRUE, showWarnings=FALSE);

#Assembly phase length variation
svg(file.path(out_folder,'controls','assembly_length_variation.svg'));
layout(matrix(c(1,3,5,2,4,6),ncol=2))
y_label = expression(paste('Assembly Rate (',min^-1,')',sep=''))
par(bty='n',mar=c(3,4.3,2,0))
data = boxplot_with_points(
	list(processed$high_Rsq$wild$results_10$assembly$slope,
		 processed$high_Rsq$S178A$results_10$assembly$slope),
	names=c('Wild-type','S178A'), notch=T, ylab=y_label, main="Minimum Phase Length: 10")
plot_dims = par("usr");
text(2.3,plot_dims[4]*0.9,paste(round(100*(data$stats[3,2]/data$stats[3,1])), '%'))

data = boxplot_with_points(
	list(processed$high_Rsq$wild$results_9$assembly$slope,
		 processed$high_Rsq$S178A$results_9$assembly$slope),
	names=c('Wild-type','S178A'), notch=T, ylab=y_label, main="Minimum Phase Length: 9")
plot_dims = par("usr");
text(2.3,plot_dims[4]*0.9,paste(round(100*(data$stats[3,2]/data$stats[3,1])), '%'))

data = boxplot_with_points(
	list(processed$high_Rsq$wild$results_8$assembly$slope,
		 processed$high_Rsq$S178A$results_8$assembly$slope),
	names=c('Wild-type','S178A'), notch=T, ylab=y_label, main="Minimum Phase Length: 8")
plot_dims = par("usr");
text(2.3,plot_dims[4]*0.9,paste(round(100*(data$stats[3,2]/data$stats[3,1])), '%'))

data = boxplot_with_points(
	list(processed$high_Rsq$wild$results_7$assembly$slope,
		 processed$high_Rsq$S178A$results_7$assembly$slope),
	names=c('Wild-type','S178A'), notch=T, ylab=y_label, main="Minimum Phase Length: 7")
plot_dims = par("usr");
text(2.3,plot_dims[4]*0.9,paste(round(100*(data$stats[3,2]/data$stats[3,1])), '%'))

data = boxplot_with_points(
	list(processed$high_Rsq$wild$results_6$assembly$slope,
		 processed$high_Rsq$S178A$results_6$assembly$slope),
	names=c('Wild-type','S178A'), notch=T, ylab=y_label, main="Minimum Phase Length: 6")
plot_dims = par("usr");
text(2.3,plot_dims[4]*0.9,paste(round(100*(data$stats[3,2]/data$stats[3,1])), '%'))

data = boxplot_with_points(
	list(processed$high_Rsq$wild$results_5$assembly$slope,
		 processed$high_Rsq$S178A$results_5$assembly$slope),
	names=c('Wild-type','S178A'), notch=T, ylab=y_label, main="Minimum Phase Length: 5")
plot_dims = par("usr");
text(2.3,plot_dims[4]*0.9,paste(round(100*(data$stats[3,2]/data$stats[3,1])), '%'))
graphics.off()

#Disassembly phase length variation
svg(file.path(out_folder,'controls','disassembly_length_variation.svg'));
layout(matrix(c(1,3,5,2,4,6),ncol=2))
y_label = expression(paste('Disassembly Rate (',min^-1,')',sep=''))
par(bty='n',mar=c(3,4.3,2,0))
data = boxplot_with_points(
	list(processed$high_Rsq$wild$results_10$disassembly$slope,
		 processed$high_Rsq$S178A$results_10$disassembly$slope),
	names=c('Wild-type','S178A'), notch=T, ylab=y_label, main="Minimum Phase Length: 10")
plot_dims = par("usr");
text(2.3,plot_dims[4]*0.9,paste(round(100*(data$stats[3,2]/data$stats[3,1])), '%'))

data = boxplot_with_points(
	list(processed$high_Rsq$wild$results_9$disassembly$slope,
		 processed$high_Rsq$S178A$results_9$disassembly$slope),
	names=c('Wild-type','S178A'), notch=T, ylab=y_label, main="Minimum Phase Length: 9")
plot_dims = par("usr");
text(2.3,plot_dims[4]*0.9,paste(round(100*(data$stats[3,2]/data$stats[3,1])), '%'))

data = boxplot_with_points(
	list(processed$high_Rsq$wild$results_8$disassembly$slope,
		 processed$high_Rsq$S178A$results_8$disassembly$slope),
	names=c('Wild-type','S178A'), notch=T, ylab=y_label, main="Minimum Phase Length: 8")
plot_dims = par("usr");
text(2.3,plot_dims[4]*0.9,paste(round(100*(data$stats[3,2]/data$stats[3,1])), '%'))

data = boxplot_with_points(
	list(processed$high_Rsq$wild$results_7$disassembly$slope,
		 processed$high_Rsq$S178A$results_7$disassembly$slope),
	names=c('Wild-type','S178A'), notch=T, ylab=y_label, main="Minimum Phase Length: 7")
plot_dims = par("usr");
text(2.3,plot_dims[4]*0.9,paste(round(100*(data$stats[3,2]/data$stats[3,1])), '%'))

data = boxplot_with_points(
	list(processed$high_Rsq$wild$results_6$disassembly$slope,
		 processed$high_Rsq$S178A$results_6$disassembly$slope),
	names=c('Wild-type','S178A'), notch=T, ylab=y_label, main="Minimum Phase Length: 6")
plot_dims = par("usr");
text(2.3,plot_dims[4]*0.9,paste(round(100*(data$stats[3,2]/data$stats[3,1])), '%'))

data = boxplot_with_points(
	list(processed$high_Rsq$wild$results_5$disassembly$slope,
		 processed$high_Rsq$S178A$results_5$disassembly$slope),
	names=c('Wild-type','S178A'), notch=T, ylab=y_label, main="Minimum Phase Length: 5")
plot_dims = par("usr");
text(2.3,plot_dims[4]*0.9,paste(round(100*(data$stats[3,2]/data$stats[3,1])), '%'))
graphics.off()


#Edge Dist at birth variation
svg(file.path(out_folder,'controls','edge_dist_at_birth.svg'));
layout(matrix(c(1,3,5,2,4,6),ncol=2))
y_label = expression(paste('Assembly Rate (',min^-1,')',sep=''))
par(bty='n',mar=c(3,4.3,2,0))
data = boxplot_with_points(
	list(processed$high_Rsq$wild$results_10$assembly$edge_dist,
		 processed$high_Rsq$S178A$results_10$assembly$edge_dist),
	names=c('Wild-type','S178A'), notch=T, ylab=y_label, main="Minimum Phase Length: 10")
plot_dims = par("usr");
text(2.3,plot_dims[4]*0.9,paste(round(100*(data$stats[3,2]/data$stats[3,1])), '%'))

data = boxplot_with_points(
	list(processed$high_Rsq$wild$results_9$assembly$edge_dist,
		 processed$high_Rsq$S178A$results_9$assembly$edge_dist),
	names=c('Wild-type','S178A'), notch=T, ylab=y_label, main="Minimum Phase Length: 9")
plot_dims = par("usr");
text(2.3,plot_dims[4]*0.9,paste(round(100*(data$stats[3,2]/data$stats[3,1])), '%'))

data = boxplot_with_points(
	list(processed$high_Rsq$wild$results_8$assembly$edge_dist,
		 processed$high_Rsq$S178A$results_8$assembly$edge_dist),
	names=c('Wild-type','S178A'), notch=T, ylab=y_label, main="Minimum Phase Length: 8")
plot_dims = par("usr");
text(2.3,plot_dims[4]*0.9,paste(round(100*(data$stats[3,2]/data$stats[3,1])), '%'))

data = boxplot_with_points(
	list(processed$high_Rsq$wild$results_7$assembly$edge_dist,
		 processed$high_Rsq$S178A$results_7$assembly$edge_dist),
	names=c('Wild-type','S178A'), notch=T, ylab=y_label, main="Minimum Phase Length: 7")
plot_dims = par("usr");
text(2.3,plot_dims[4]*0.9,paste(round(100*(data$stats[3,2]/data$stats[3,1])), '%'))

data = boxplot_with_points(
	list(processed$high_Rsq$wild$results_6$assembly$edge_dist,
		 processed$high_Rsq$S178A$results_6$assembly$edge_dist),
	names=c('Wild-type','S178A'), notch=T, ylab=y_label, main="Minimum Phase Length: 6")
plot_dims = par("usr");
text(2.3,plot_dims[4]*0.9,paste(round(100*(data$stats[3,2]/data$stats[3,1])), '%'))

data = boxplot_with_points(
	list(processed$high_Rsq$wild$results_5$assembly$edge_dist,
		 processed$high_Rsq$S178A$results_5$assembly$edge_dist),
	names=c('Wild-type','S178A'), notch=T, ylab=y_label, main="Minimum Phase Length: 5")
plot_dims = par("usr");
text(2.3,plot_dims[4]*0.9,paste(round(100*(data$stats[3,2]/data$stats[3,1])), '%'))
graphics.off()


#Edge dist at death variation
svg(file.path(out_folder,'controls','edge_dist_at_death.svg'));
layout(matrix(c(1,3,5,2,4,6),ncol=2))
y_label = expression(paste('Disassembly Rate (',min^-1,')',sep=''))
par(bty='n',mar=c(3,4.3,2,0))
data = boxplot_with_points(
	list(processed$high_Rsq$wild$results_10$disassembly$edge_dist,
		 processed$high_Rsq$S178A$results_10$disassembly$edge_dist),
	names=c('Wild-type','S178A'), notch=T, ylab=y_label, main="Minimum Phase Length: 10")
plot_dims = par("usr");
text(2.3,plot_dims[4]*0.9,paste(round(100*(data$stats[3,2]/data$stats[3,1])), '%'))

data = boxplot_with_points(
	list(processed$high_Rsq$wild$results_9$disassembly$edge_dist,
		 processed$high_Rsq$S178A$results_9$disassembly$edge_dist),
	names=c('Wild-type','S178A'), notch=T, ylab=y_label, main="Minimum Phase Length: 9")
plot_dims = par("usr");
text(2.3,plot_dims[4]*0.9,paste(round(100*(data$stats[3,2]/data$stats[3,1])), '%'))

data = boxplot_with_points(
	list(processed$high_Rsq$wild$results_8$disassembly$edge_dist,
		 processed$high_Rsq$S178A$results_8$disassembly$edge_dist),
	names=c('Wild-type','S178A'), notch=T, ylab=y_label, main="Minimum Phase Length: 8")
plot_dims = par("usr");
text(2.3,plot_dims[4]*0.9,paste(round(100*(data$stats[3,2]/data$stats[3,1])), '%'))

data = boxplot_with_points(
	list(processed$high_Rsq$wild$results_7$disassembly$edge_dist,
		 processed$high_Rsq$S178A$results_7$disassembly$edge_dist),
	names=c('Wild-type','S178A'), notch=T, ylab=y_label, main="Minimum Phase Length: 7")
plot_dims = par("usr");
text(2.3,plot_dims[4]*0.9,paste(round(100*(data$stats[3,2]/data$stats[3,1])), '%'))

data = boxplot_with_points(
	list(processed$high_Rsq$wild$results_6$disassembly$edge_dist,
		 processed$high_Rsq$S178A$results_6$disassembly$edge_dist),
	names=c('Wild-type','S178A'), notch=T, ylab=y_label, main="Minimum Phase Length: 6")
plot_dims = par("usr");
text(2.3,plot_dims[4]*0.9,paste(round(100*(data$stats[3,2]/data$stats[3,1])), '%'))

data = boxplot_with_points(
	list(processed$high_Rsq$wild$results_5$disassembly$edge_dist,
		 processed$high_Rsq$S178A$results_5$disassembly$edge_dist),
	names=c('Wild-type','S178A'), notch=T, ylab=y_label, main="Minimum Phase Length: 5")
plot_dims = par("usr");
text(2.3,plot_dims[4]*0.9,paste(round(100*(data$stats[3,2]/data$stats[3,1])), '%'))
graphics.off()
