library(Biobase)
library(data.table)
library(pcaMethods)
library(doParallel)
suppressWarnings(library(EnDecon))

set("./simulation/run methods/")

### python env
python_env <- "~/.conda/envs/EnDecon_GPU/bin/python"
sim_pbac_seg_muraro_wang <- readRDS("./simulation/data/sim_pbac_seg_muraro_wang.Rds")

######################################## non-batch(Scenario 1) ########################################
Results.dec.ori <- list()
sc_ref_ori <- as.matrix(sim_pbac_seg_muraro_wang$ref_data$sc_exp)
sc_ref_label_ori <- sim_pbac_seg_muraro_wang$ref_data$sc_label
names(sc_ref_label_ori) <- colnames(sc_ref_ori)
spot_loc_seg <- as.data.frame(sim_pbac_seg_muraro_wang$spot_loc)
for (i in 1:length(sim_pbac_seg_muraro_wang$spot_data)) {
  cat("i th:", i, "\n")

  temp <- sim_pbac_seg_muraro_wang$spot_data[[i]]
  spot_exp <- as.matrix(temp$sim_spots)
  colnames(spot_exp) <- rownames(spot_loc_seg)

  Results.dec.ori[[i]] <- EnDecon_individual_methods(sc_exp = sc_ref_ori, sc_label = sc_ref_label_ori,
                                                     spot_exp = spot_exp, spot_loc = spot_loc_seg,
                                                     python_env = python_env, use_gpu = TRUE,
                                                     gene_det_in_min_cells_per = 0.01,
                                                     RCTD.CELL_MIN_INSTANCE = 5, saving_results = FALSE)
}
# save(Results.dec.ori, file = "Results.dec.ori.RData")
### ensemble
Results.dec.ori.decon <- list()
# Results.dec.muraro.times <- list()
Results.decon.ori.mean <- list()
Results.dec.ori.weighte <- list()
list_w_ensemble_ori <- list()


for (i in 1:length(Results.dec.ori)) {
  cat("i th: ", i ,"\n")

  Results.Deconv <- Results.dec.ori[[i]][[1]]
  Results.times <- Results.dec.ori[[i]][[2]]
  ### weight ensemble
  pred.weight.ensemble <- Sys.time()
  Results <- solve_ensemble(Results.Deconv = Results.Deconv, lambda = NULL, prob.quantile = 0.5)
  end.weight.ensemble <- Sys.time()
  time.weight.ensemble <- difftime(end.weight.ensemble, pred.weight.ensemble, units = "mins")
  Results.dec.ori.weighte[[i]] <- Results$V_norm
  list_w_ensemble_ori[[i]] <- Results$w
  ### mean ensemble
  pred.mean.ensemble <- Sys.time()
  Results.decon.ori.mean[[i]] <- Reduce("+", Map("*", Results.Deconv, 1))/length(Results.Deconv)
  end.mean.ensemble <- Sys.time()
  time.mean.ensemble <- difftime(end.mean.ensemble, pred.mean.ensemble, units = "mins")
  ### ensemble setting
  Results.dec.ori[[i]][[1]]$EnDecon_mean <- Results.decon.ori.mean[[i]]
  Results.dec.ori[[i]][[1]]$EnDecon <- Results$V_norm
  ### ensemble times
  Results.dec.ori[[i]][[2]] <- c(Results.dec.ori[[i]][[2]], time.mean.ensemble)
  Results.dec.ori[[i]][[2]] <- c(Results.dec.ori[[i]][[2]], time.weight.ensemble)
  ### resort the name of methods
  Results.dec.ori[[i]][[1]] <- Results.dec.ori[[i]][[1]][method_index]
  names(Results.dec.ori[[i]][[1]]) <- method_names
  Results.dec.ori[[i]][[2]] <- Results.dec.ori[[i]][[2]][method_index]
  names(Results.dec.ori[[i]][[2]]) <- method_names

  ### resort the name of weight
  list_w_ensemble_ori[[i]] <- list_w_ensemble_ori[[i]][method_index[1:14]]
  names(list_w_ensemble_ori[[i]]) <- method_names[1:14]
}
save(Results.dec.ori, file = "Results.dec.ori.all.RData")

################################## Segerstolpe reference (Scenario 2) ####################################
Results.dec.seg <- list()
sc_ref_seg <- as.matrix(sim_pbac_seg_muraro_wang$ref_batch_seg$sc_exp)
sc_ref_label_seg <- sim_pbac_seg_muraro_wang$ref_batch_seg$sc_label
names(sc_ref_label_seg) <- colnames(sc_ref_seg)
spot_loc_seg <- as.data.frame(sim_pbac_seg_muraro_wang$spot_loc)
for (i in 1:length(sim_pbac_seg_muraro_wang$spot_data)) {
  cat("i th:", i, "\n")
  temp <- sim_pbac_seg_muraro_wang$spot_data[[i]]
  spot_exp <- as.matrix(temp$sim_spots)
  colnames(spot_exp) <- rownames(spot_loc_seg)
  Results.dec.seg[[i]] <- EnDecon_individual_methods(sc_exp = sc_ref_seg, sc_label = sc_ref_label_seg,
                                                     spot_exp = spot_exp, spot_loc = spot_loc_seg,
                                                     python_env = python_env, use_gpu = TRUE,
                                                     gene_det_in_min_cells_per = 0.01,
                                                     RCTD.CELL_MIN_INSTANCE = 5, saving_results = FALSE)
}
# save(Results.dec.seg, file = "Results.dec.seg.RData")
### ensemble
Results.dec.seg.decon <- list()
Results.decon.seg.mean <- list()
Results.dec.seg.weighte <- list()
list_w_ensemble_seg <- list()
for (i in 1:length(Results.dec.seg)) {
  cat("i th: ", i ,"\n")
  Results.Deconv <- Results.dec.seg[[i]][[1]]
  Results.times <- Results.dec.seg[[i]][[2]]
  ### weight ensemble
  pred.weight.ensemble <- Sys.time()
  Results <- solve_ensemble(Results.Deconv = Results.Deconv, lambda = NULL, prob.quantile = 0.5)
  end.weight.ensemble <- Sys.time()
  time.weight.ensemble <- difftime(end.weight.ensemble, pred.weight.ensemble, units = "mins")
  Results.dec.seg.weighte[[i]] <- Results$V_norm
  list_w_ensemble_seg[[i]] <- Results$w
  ### mean ensemble
  pred.mean.ensemble <- Sys.time()
  Results.decon.seg.mean[[i]] <- Reduce("+", Map("*", Results.Deconv, 1))/length(Results.Deconv)
  end.mean.ensemble <- Sys.time()
  time.mean.ensemble <- difftime(end.mean.ensemble, pred.mean.ensemble, units = "mins")
  ### ensemble setting
  Results.dec.seg[[i]][[1]]$EnDecon_mean <- Results.decon.seg.mean[[i]]
  Results.dec.seg[[i]][[1]]$EnDecon <- Results$V_norm
  ### ensemble times
  Results.dec.seg[[i]][[2]] <- c(Results.dec.seg[[i]][[2]], time.mean.ensemble)
  Results.dec.seg[[i]][[2]] <- c(Results.dec.seg[[i]][[2]], time.weight.ensemble)
  ### resort the name of methods
  Results.dec.seg[[i]][[1]] <- Results.dec.seg[[i]][[1]][method_index]
  names(Results.dec.seg[[i]][[1]]) <- method_names
  Results.dec.seg[[i]][[2]] <- Results.dec.seg[[i]][[2]][method_index]
  names(Results.dec.seg[[i]][[2]]) <- method_names
  ### add weight
  list_w_ensemble_seg[[i]] <- list_w_ensemble_seg[[i]][method_index[1:14]]
  names(list_w_ensemble_seg[[i]]) <- method_names[1:14]

  list_w_ensemble_seg[[i]] <- list_w_ensemble_seg[[i]][method_index[1:14]]
  names(list_w_ensemble_seg[[i]]) <- method_names[1:14]
}
save(Results.dec.seg, file = "Results.dec.seg.all.RData")
################################## Segerstolpe reference (Scenario 2) ####################################

########################################### Muraro reference (Scenario 2)  ###############################
Results.dec.muraro <- list()
sc_ref_muraro_temp <- t(as.matrix(sim_pbac_seg_muraro_wang$ref_batch_muraro$sc_exp))
sc_ref_muraro <- matrix(as.integer(sc_ref_muraro_temp), nrow = nrow(sc_ref_muraro_temp), ncol = ncol(sc_ref_muraro_temp))
rownames(sc_ref_muraro) <- rownames(sc_ref_muraro_temp)
colnames(sc_ref_muraro) <- colnames(sc_ref_muraro_temp)
sc_ref_label_muraro <- sim_pbac_seg_muraro_wang$ref_batch_muraro$sc_label
names(sc_ref_label_muraro) <- colnames(sc_ref_muraro)
spot_loc_seg <- as.data.frame(sim_pbac_seg_muraro_wang$spot_loc)
for (i in 1:length(sim_pbac_seg_muraro_wang$spot_data)) {
  cat("i th:", i, "\n")
  temp <- sim_pbac_seg_muraro_wang$spot_data[[i]]
  spot_exp <- as.matrix(temp$sim_spots)
  colnames(spot_exp) <- rownames(spot_loc_seg)

  Results.dec.muraro[[i]] <- EnDecon_individual_methods(sc_exp = sc_ref_muraro, sc_label = sc_ref_label_muraro,
                                                        spot_exp = spot_exp, spot_loc = spot_loc_seg,
                                                        python_env = python_env, use_gpu = TRUE,
                                                        gene_det_in_min_cells_per = 0.01,
                                                        RCTD.CELL_MIN_INSTANCE = 5, saving_results = FALSE)
}
# save(Results.dec.muraro, file = "Results.dec.muraro.RData")
Results.dec.muraro.decon <- list()
# Results.dec.muraro.times <- list()
Results.decon.muraro.mean <- list()
Results.dec.muraro.weighte <- list()
list_w_ensemble_muraro <- list()
for (i in 1:length(Results.dec.muraro)) {
  cat("i th: ", i ,"\n")

  Results.Deconv <- Results.dec.muraro[[i]][[1]]
  Results.times <- Results.dec.muraro[[i]][[2]]
  ### weight ensemble
  pred.weight.ensemble <- Sys.time()
  Results <- solve_ensemble_norm(Results.Deconv = Results.Deconv, lambda = NULL, prob.quantile = 0.5)
  end.weight.ensemble <- Sys.time()
  time.weight.ensemble <- difftime(end.weight.ensemble, pred.weight.ensemble, units = "mins")
  Results.dec.muraro.weighte[[i]] <- Results$V_norm
  list_w_ensemble_muraro[[i]] <- Results$w
  ### mean ensemble
  pred.mean.ensemble <- Sys.time()
  Results.decon.muraro.mean[[i]] <- Reduce("+", Map("*", Results.Deconv, 1))/length(Results.Deconv)
  end.mean.ensemble <- Sys.time()
  time.mean.ensemble <- difftime(end.mean.ensemble, pred.mean.ensemble, units = "mins")
  ### ensemble setting
  Results.dec.muraro[[i]][[1]]$EnDecon_mean <- Results.decon.muraro.mean[[i]]
  Results.dec.muraro[[i]][[1]]$EnDecon <- Results$V_norm
  ### ensemble times
  Results.dec.muraro[[i]][[2]] <- c(Results.dec.muraro[[i]][[2]], time.mean.ensemble)
  Results.dec.muraro[[i]][[2]] <- c(Results.dec.muraro[[i]][[2]], time.weight.ensemble)
  ### resort the name of methods
  Results.dec.muraro[[i]][[1]] <- Results.dec.muraro[[i]][[1]][method_index]
  names(Results.dec.muraro[[i]][[1]]) <- method_names
  Results.dec.muraro[[i]][[2]] <- Results.dec.muraro[[i]][[2]][method_index]
  names(Results.dec.muraro[[i]][[2]]) <- method_names

  ### resort the name of weight
  list_w_ensemble_muraro[[i]] <- list_w_ensemble_muraro[[i]][method_index[1:14]]
  names(list_w_ensemble_muraro[[i]]) <- method_names[1:14]
}
save(Results.dec.muraro, file = "Results.dec.muraro.all.RData")
########################################### Muraro reference (Scenario 2)  ###############################
########################################### Wang reference (Scenario 2) ####################################
Results.dec.wang <- list()
sc_ref_wang_temp <- t(as.matrix(sim_pbac_seg_muraro_wang$ref_batch_wang$sc_exp))
sc_ref_wang <- matrix(as.integer(sc_ref_wang_temp), nrow = nrow(sc_ref_wang_temp), ncol = ncol(sc_ref_wang_temp))
rownames(sc_ref_wang) <- rownames(sc_ref_wang_temp)
colnames(sc_ref_wang) <- colnames(sc_ref_wang_temp)
sc_ref_wang <- t(sc_ref_wang)
sc_ref_label_wang <- sim_pbac_seg_muraro_wang$ref_batch_wang$sc_label
names(sc_ref_label_wang) <- colnames(sc_ref_wang)
spot_loc_seg <- as.data.frame(sim_pbac_seg_muraro_wang$spot_loc)

for (i in 1:length(sim_pbac_seg_muraro_wang$spot_data)) {
  cat("i th:", i, "\n")

  temp <- sim_pbac_seg_muraro_wang$spot_data[[i]]
  spot_exp <- as.matrix(temp$sim_spots)
  colnames(spot_exp) <- rownames(spot_loc_seg)

  Results.dec.wang[[i]] <- EnDecon_individual_methods(sc_exp = sc_ref_wang, sc_label = sc_ref_label_wang,
                                                      spot_exp = spot_exp, spot_loc = spot_loc_seg,
                                                      python_env = python_env, use_gpu = TRUE,
                                                      gene_det_in_min_cells_per = 0.01,
                                                      RCTD.CELL_MIN_INSTANCE = 5, saving_results = FALSE)
}
# save(Results.dec.wang, file = "Results.dec.wang.RData")
Results.dec.wang.decon <- list()
# Results.dec.muraro.times <- list()
Results.decon.wang.mean <- list()
Results.dec.wang.weighte <- list()
list_w_ensemble_wang <- list()
for (i in 1:length(Results.dec.wang)) {
  cat("i th: ", i ,"\n")

  Results.Deconv <- Results.dec.wang[[i]][[1]]
  Results.times <- Results.dec.wang[[i]][[2]]
  ### weight ensemble
  pred.weight.ensemble <- Sys.time()
  Results <- solve_ensemble_norm(Results.Deconv = Results.Deconv, lambda = NULL, prob.quantile = 0.5)
  end.weight.ensemble <- Sys.time()
  time.weight.ensemble <- difftime(end.weight.ensemble, pred.weight.ensemble, units = "mins")
  Results.dec.wang.weighte[[i]] <- Results$V_norm
  list_w_ensemble_wang[[i]] <- Results$w
  ### mean ensemble
  pred.mean.ensemble <- Sys.time()
  Results.decon.wang.mean[[i]] <- Reduce("+", Map("*", Results.Deconv, 1))/length(Results.Deconv)
  end.mean.ensemble <- Sys.time()
  time.mean.ensemble <- difftime(end.mean.ensemble, pred.mean.ensemble, units = "mins")
  ### ensemble setting
  Results.dec.wang[[i]][[1]]$EnDecon_mean <- Results.decon.wang.mean[[i]]
  Results.dec.wang[[i]][[1]]$EnDecon <- Results$V_norm
  ### ensemble times
  Results.dec.wang[[i]][[2]] <- c(Results.dec.wang[[i]][[2]], time.mean.ensemble)
  Results.dec.wang[[i]][[2]] <- c(Results.dec.wang[[i]][[2]], time.weight.ensemble)
  ### resort the name of methods
  Results.dec.wang[[i]][[1]] <- Results.dec.wang[[i]][[1]][method_index]
  names(Results.dec.wang[[i]][[1]]) <- method_names
  Results.dec.wang[[i]][[2]] <- Results.dec.wang[[i]][[2]][method_index]
  names(Results.dec.wang[[i]][[2]]) <- method_names

  list_w_ensemble_wang[[i]] <- list_w_ensemble_wang[[i]][method_index[1:14]]
  names(list_w_ensemble_wang[[i]]) <- method_names[1:14]
}
save(Results.dec.wang, file = "Results.dec.wang.all.RData")
########################################### Wang reference (Scenario 2) ####################################


