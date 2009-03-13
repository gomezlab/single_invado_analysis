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

out_folder = '../../doc/publication/figures/'

library(boot)

files_to_load = c('../single_lin.csv', '../edge_velo_pre_birth.csv', 
				  '../edge_velo_post_birth.csv', '../edge_velo_null.csv', 
				  '../edge_velo_post_birth.csv', '../edge_velo_pre_death.csv')

raw <- load_data_files(exp_dirs, files_to_load, c(T,F,F,F, F, F));
raw <- lapply(raw, function(x) x[2:10])

lineage_data <- raw[[1]]
pre_birth_raw <- raw[[2]]
post_death_raw <- raw[[3]]
null_raw <- raw[[4]]
post_birth_raw <- raw[[5]]
pre_death_raw <- raw[[6]]

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

	rows_to_include = apply(this_pre_birth, 1, function(x) all(x[(length(x) - 2):length(x)] < 0))
	rows_to_include[is.na(rows_to_include)] = FALSE
	pre_birth$protrude_means[[i]] = as.numeric(mean(this_pre_birth[rows_to_include,], na.rm=T))
	
	temp = find_col_conf_ints(this_pre_birth[rows_to_include,])
	pre_birth$protrude_upper[[i]] = temp$u
	pre_birth$protrude_lower[[i]] = temp$l

}

nums = c()
pos_count = c()
neg_count = c()
exp_num = 1
for (i in 0:55) {
	this_data = as.matrix(pre_birth$data[[exp_num]][,length(pre_birth_raw[[exp_num]]) - i])
	nums = c(nums, sum(! is.na(this_data)))

	mat_data = as.matrix(pre_birth$data[[exp_num]][,(length(pre_birth_raw[[exp_num]]) - i):length(pre_birth_raw[[exp_num]])])
	all_pos = apply(mat_data, 1, function(x) all(x > 0))
	all_pos[is.na(all_pos)] = FALSE
	pos_count = c(pos_count, sum(all_pos))

	all_neg = apply(mat_data, 1, function(x) all(x < 0))
	all_neg[is.na(all_neg)] = FALSE
	neg_count = c(neg_count, sum(all_neg))
}
pos_count = pos_count/nums
neg_count = neg_count/nums

pos_null = c()
neg_null = c()
for (i in 1:56) {
	pos_percent = sum(pre_birth$unlist[[exp_num]] > 0)/length(pre_birth$unlist[[exp_num]])
	neg_percent = sum(pre_birth$unlist[[exp_num]] < 0)/length(pre_birth$unlist[[exp_num]])
	
	pos_null = c(pos_null, (pos_percent^i*nums[i])/nums[i])
	neg_null = c(neg_null, (neg_percent^i*nums[i])/nums[i])
}

pdf(file.path(out_folder,'edge_experiment_plots.pdf'))
par(mfcol=c(3,3),bty='n',mar=c(4.2,4.1,2,0.2), oma=c(0,0,2,0))
for (i in 1:length(pre_birth$means)) {
	y_range = c()
	y_range[1] = min(c(pre_birth$lower[[i]], pre_birth$protrude_lower[[i]]), na.rm=T);
	y_range[2] = max(c(pre_birth$upper[[i]], pre_birth$protrude_upper[[i]]), na.rm=T);
	errbar(1:length(pre_birth$means[[i]]), pre_birth$means[[i]], pre_birth$upper[[i]], pre_birth$lower[[i]], ylab = i, xlab = "", ylim=y_range)

	segments(1, pre_birth$all_u[[i]], length(pre_birth$means[[i]]), pre_birth$all_u[[i]], col='green')
	segments(1, pre_birth$all_l[[i]], length(pre_birth$means[[i]]), pre_birth$all_l[[i]], col='green')
	
	errbar(1:length(pre_birth$means[[i]]), pre_birth$protrude_means[[i]], pre_birth$protrude_upper[[i]], pre_birth$protrude_lower[[i]], col='red', add=T)

	if (i == 1) {	
		mtext('Pre-birth', outer=T)
	}
}

for (i in 1:length(post_birth$means)) {
	errbar(1:length(post_birth$means[[i]]), post_birth$means[[i]], post_birth$upper[[i]], post_birth$lower[[i]], ylab = i, xlab = "")

	segments(1, post_birth$all_u[[i]], length(post_birth$means[[i]]), post_birth$all_u[[i]], col='green')
	segments(1, post_birth$all_l[[i]], length(post_birth$means[[i]]), post_birth$all_l[[i]], col='green')
	if (i == 1) {
		mtext('Post-birth', outer=T)
	}
}

for (i in 1:length(pre_death$means)) {
	errbar(1:length(pre_death$means[[i]]), pre_death$means[[i]], pre_death$upper[[i]], pre_death$lower[[i]], ylab = i, xlab = "")

	segments(1, pre_death$all_u[[i]], length(pre_birth$means[[i]]), pre_death$all_u[[i]], col='green')
	segments(1, pre_death$all_l[[i]], length(pre_birth$means[[i]]), pre_death$all_l[[i]], col='green')

	if (i == 1) {
		mtext('Pre-death', outer=T)
	}
}

for (i in 1:length(post_death$means)) {
	errbar(1:length(post_death$means[[i]]), post_death$means[[i]], post_death$upper[[i]], post_death$lower[[i]], ylab = i, xlab = "")
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
