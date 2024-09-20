#devtools::install_github("rcannood/SCORPIUS")
library(SCORPIUS)

endo<- readRDS("result/endothelial_cells_rmbatch.rds")
####至此，github官网上expression用的是data的数据，而不是scale data，试试看有什么区别
expression_data <- t(as.matrix(endo@assays$RNA@data))
expression_scale <- t(as.matrix(endo@assays$RNA@scale.data))

group_name =  as.factor(as.character(endo$group))
table(group_name)
dim(expression_data)
dim(expression_scale)

# #######expression_data----
# space <- reduce_dimensionality(expression_data, "spearman")
# draw_trajectory_plot(space, group_name, contour = TRUE,)
# traj <- infer_trajectory(space)
# draw_trajectory_plot(space, group_name, traj$path, contour = TRUE)
# # warning: setting num_permutations to 10 requires a long time (~30min) to run!
# # set it to 0 and define a manual cutoff for the genes (e.g. top 200) for a much shorter execution time.
# gimp <- gene_importances(
#   expression_data,
#   traj$time,
#   num_permutations = 0,
#   num_threads = 8,
#   ntree = 10000,
#   ntree_perm = 1000
# )
# #gimp$qvalue <- p.adjust(gimp$pvalue, "BH", length(gimp$pvalue))
# gene_sel <- gimp$gene[1:200]
# expr_sel <- scale_quantile(expression_data[,gene_sel])
# 
# # Draw a time series heatmap
# time <- traj$time
# draw_trajectory_heatmap(expr_sel, time)
# 
# ## Also show the progression groupings
# draw_trajectory_heatmap(expr_sel, time,
#                         progression_group=group_name)

######expression_scale----
set.seed(111)
space2 <- reduce_dimensionality(expression_scale, "spearman")
draw_trajectory_plot(space2, group_name, contour = TRUE,)
traj2 <- infer_trajectory(space2)
draw_trajectory_plot(space2, group_name, traj2$path, contour = TRUE)
# warning: setting num_permutations to 10 requires a long time (~30min) to run!
# set it to 0 and define a manual cutoff for the genes (e.g. top 200) for a much shorter execution time.
gimp2 <- gene_importances(
  expression_scale, 
  traj2$time, 
  num_permutations = 0, 
  num_threads = 8, 
  ntree = 10000,
  ntree_perm = 1000
) 
#gimp2$qvalue <- p.adjust(gimp2$pvalue, "BH", length(gimp2$pvalue))
gene_sel2 <- gimp2$gene[1:200]
expr_sel2 <- scale_quantile(expression_scale[,gene_sel2])

# Draw a time series heatmap
time2 <- traj2$time
# draw_trajectory_heatmap(expr_sel2, time2)
# 
# ## Also show the progression groupings
# draw_trajectory_heatmap(expr_sel2, time2, 
#                         progression_group=group_name)
# draw_trajectory_heatmap(expr_sel2, time2, 
#                         progression_group=group_name,
#                         show_labels_row=T)
modules_2 <- extract_modules(scale_quantile(expr_sel2), traj2$time, verbose = F)

pdf("plot_new/SCORPIUS_endo.pdf")
draw_trajectory_plot(space2, group_name, traj2$path, contour = TRUE)
draw_trajectory_heatmap(expr_sel2, time2, 
                        progression_group=group_name,
                        show_labels_row=T, fontsize_row = 4,
                        modules=modules_2)
dev.off()
