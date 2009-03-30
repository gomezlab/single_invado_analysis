rm(list = ls())
source('FA_analysis_lib.R')
library(lattice)
library(geneplotter)
library(Hmisc)

################################################################################
#Result loading
################################################################################
exp_dirs <- Sys.glob('../../results/focal_adhesions/*/adhesion_props/lin_time_series/')
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]

out_folder = '../../doc/publication/figures/'

library(boot)

files_to_load = c('../single_lin.csv', '../Edge_projection/pre_birth.csv', 
				  '../Edge_projection/post_birth.csv', '../Edge_projection/null.csv', 
				  '../Edge_projection/pre_death.csv', '../Edge_projection/post_death.csv')

raw <- load_data_files(exp_dirs, files_to_load, c(T, F, F, F, F, F));
raw <- lapply(raw, function(x) x[2:10])

lineage_data <- raw[[1]]
pre_birth_raw <- raw[[2]]
post_birth_raw <- raw[[3]]
null_raw <- raw[[4]]
pre_death_raw <- raw[[5]]
post_death_raw <- raw[[6]]

exp_names <- raw[[length(raw)]]

########################################
#Velocity Rings
########################################
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


################################################################################
# Edge Velocity Processing
################################################################################

pre_birth = list()
post_birth = list()
pre_death = list()
post_death = list()
null = list()

all_pre_birth = matrix(NA, nrow = 1, ncol = 20)
for (i in 1:length(pre_birth_raw)) {
	this_pre_birth = pre_birth_raw[[i]];
	this_post_death = post_death_raw[[i]];
	this_null_raw = null_raw[[i]];
	this_post_birth = post_birth_raw[[i]];
	this_pre_death = pre_death_raw[[i]];
	
	birth_filt = ! lineage_data[[i]]$split_birth_status & 
				 ! is.na(this_pre_birth[,length(this_pre_birth)]) & 
				 lineage_data[[i]]$starting_edge_dist < Inf
	print(sum(birth_filt)/length(lineage_data[[i]]$starting_edge_dist))

	death_filt = as.logical(lineage_data[[i]]$death_status) 
	
	this_pre_birth = this_pre_birth[birth_filt,]
	pre_birth$data[[i]] = this_pre_birth
	this_post_birth = this_post_birth[birth_filt,]
	post_birth$data[[i]] = this_post_birth
	this_pre_death = this_pre_death[death_filt,] 
	pre_death$data[[i]] = this_pre_death
	this_post_death = this_post_death[death_filt,] 
	post_death$data[[i]] = this_post_death
	
	all_pre_birth = rbind(all_pre_birth, as.matrix(this_pre_birth[,(dim(this_pre_birth)[[2]] - 19):dim(this_pre_birth)[[2]]]))
	
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

	rows_to_include = apply(this_pre_birth, 1, function(x) all(x[(length(x) - 2):length(x)] < 0))
	rows_to_include[is.na(rows_to_include)] = FALSE
	
	temp = find_col_conf_ints(this_pre_birth[rows_to_include,])
	pre_birth$protrude_upper[[i]] = temp$u
	pre_birth$protrude_lower[[i]] = temp$l
	pre_birth$protrude_means[[i]] = temp$mean_vals

	temp = find_col_conf_ints(this_post_birth[rows_to_include,])
	post_birth$protrude_upper[[i]] = temp$u
	post_birth$protrude_lower[[i]] = temp$l
	post_birth$protrude_means[[i]] = temp$mean_vals
}

################################################################################
# Plotting
################################################################################

pdf(file.path(out_folder,'edge_experiment_plots.pdf'))
par(mfcol=c(3,3),bty='n',mar=c(4.2,4.1,2,0.2), oma=c(0,0,2,0))
for (i in 1:length(pre_birth$means)) {
	y_range = c()
	y_range[1] = min(c(pre_birth$lower[[i]], pre_birth$protrude_lower[[i]]), na.rm=T);
	y_range[2] = max(c(pre_birth$upper[[i]], pre_birth$protrude_upper[[i]]), na.rm=T);
	
	x_range = (length(pre_birth$means[[i]])- 9):length(pre_birth$means[[i]]);
	
	errbar(x_range, 
		   pre_birth$means[[i]][x_range], 
		   pre_birth$upper[[i]][x_range], 
		   pre_birth$lower[[i]][x_range], 
		   main = exp_names[[i]], ylab = "", xlab = "")

#	segments(1, pre_birth$all_u[[i]], length(pre_birth$means[[i]]), pre_birth$all_u[[i]], col='green')
#	segments(1, pre_birth$all_l[[i]], length(pre_birth$means[[i]]), pre_birth$all_l[[i]], col='green')
#	
#	errbar(1:length(pre_birth$means[[i]]), pre_birth$protrude_means[[i]], pre_birth$protrude_upper[[i]], pre_birth$protrude_lower[[i]], col='red', add=T)

	if (i == 1) {	
		mtext('Pre-birth', outer=T)
	}
}

for (i in 1:length(post_birth$means)) {
	y_range = c()
	y_range[1] = min(c(post_birth$lower[[i]], post_birth$protrude_lower[[i]]), na.rm=T);
	y_range[2] = max(c(post_birth$upper[[i]], post_birth$protrude_upper[[i]]), na.rm=T);
	
	x_range = 1:10;
	
	errbar(x_range, 
		   post_birth$means[[i]][x_range], 
		   post_birth$upper[[i]][x_range], 
		   post_birth$lower[[i]][x_range], 
		   ylab = i, xlab = "", ylim=y_range)

	segments(1, post_birth$all_u[[i]], length(post_birth$means[[i]]), post_birth$all_u[[i]], col='green')
	segments(1, post_birth$all_l[[i]], length(post_birth$means[[i]]), post_birth$all_l[[i]], col='green')
	
	errbar(1:length(post_birth$means[[i]]), 
		   post_birth$protrude_means[[i]], 
		   post_birth$protrude_upper[[i]], 
		   post_birth$protrude_lower[[i]], 
		   col='red', add=T)

	
	if (i == 1) {
		mtext('Post-birth', outer=T)
	}
}

for (i in 1:length(pre_death$means)) {
	
	x_range = (length(pre_death$means[[i]])- 9):length(pre_death$means[[i]]);

	errbar(x_range, 
		   pre_death$means[[i]][x_range], 
		   pre_death$upper[[i]][x_range], 
		   pre_death$lower[[i]][x_range], 
		   ylab = i, xlab = "")

	segments(1, pre_death$all_u[[i]], length(pre_birth$means[[i]]), pre_death$all_u[[i]], col='green')
	segments(1, pre_death$all_l[[i]], length(pre_birth$means[[i]]), pre_death$all_l[[i]], col='green')

	if (i == 1) {
		mtext('Pre-death', outer=T)
	}
}

for (i in 1:length(post_death$means)) {
	
	x_range = 1:10;
	
	errbar(x_range, 
		   post_death$means[[i]][x_range], 
		   post_death$upper[[i]][x_range], 
		   post_death$lower[[i]][x_range], 
		   ylab = i, xlab = "")
		   
	segments(1, post_death$all_u[[i]], length(pre_birth$means[[i]]), post_death$all_u[[i]], col='green')
	segments(1, post_death$all_l[[i]], length(pre_birth$means[[i]]), post_death$all_l[[i]], col='green')
	
	if (i == 1) {
		mtext('Post-death', outer=T)
	}
}

for (i in 1:length(pre_birth$unlist)) {
	hist(pre_birth$unlist[[i]], main=i,xlab='Projected Edge Velocity')
	if (i == 1) {
		mtext('Pre-birth Data Sets', outer=T)
	}	
}

graphics.off()