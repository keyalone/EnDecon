library(ggplot2)
sim_pbac_seg_muraro_wang <- readRDS("./simulation/data/sim_pbac_seg_muraro_wang.Rds")
method_index <- c(1, 2, 3, 4, 5, 8, 7,  9, 10, 11, 12, 14, 13, 6, 15, 16)
method_names <- c("CARD", "Cell2location", "DeconRNASeq", "DestVI", "DWLS", "MuSiC all gene",
                  "MuSiC weighted", "RCTD", "SCDC", "SpatialDWLS", "SPOTlight",
                  "STdeconvolve", "Stereoscope", "SVR", "EnDecon_mean", "EnDecon")
load("./simulation/run methods/Results.dec.seg.all.RData")
plot.list.seg <- list()
library(philentropy)
library(cowplot)
niter <- 10
R_square_seg = matrix(NA, nrow = niter, ncol = length(method_index))
RMSE_seg = matrix(NA,     nrow =niter, ncol = length(method_index))
cell_type_seg <- sim_pbac_seg_muraro_wang$ref_batch_seg$cell_type
for (k in 1:length(method_index)) {
  for (i in 1:niter) {
    temp1 <- Results.dec.seg[[i]][[1]][[k]][,cell_type_seg]
    temp2 <- sim_pbac_seg_muraro_wang$spot_data[[i]]$sim_spots_pro[,cell_type_seg]
    R_square_seg[i, k] <- cor(c(temp1), c(temp2), method = "pearson")
    RMSE_seg[i, k] <- sqrt(mean((temp2- temp1)^2))
  }
}
colnames(R_square_seg) = method_names
colnames(RMSE_seg) = method_names
sort(apply(R_square_seg,2, median))
sort(apply(RMSE_seg,2, median))

####barplot in R

df1_box <- data.frame(PCC = c(R_square_seg),
                      RMSE = c(RMSE_seg),
                      Method = factor(rep(method_names,each=niter),
                                      levels = method_names))

error_box = ggplot(df1_box, aes(x = Method , y = PCC,
                                fill = Method))+
  geom_boxplot()+
  theme_classic()+
  labs(y="PCC",x=" ", title = "Scenario 2")+
  theme(
    axis.text.x.bottom = element_blank(),
    axis.text.y.left = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.ticks.x = element_blank(),
    plot.title = element_text( size=14,hjust = 0.5),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")
  )
# error_box
plot.list.seg[[1]] <- error_box

RMSE_box = ggplot(df1_box, aes(x = Method , y = RMSE,
                               fill = Method))+
  geom_boxplot()+
  theme_classic()+
  labs(y="RMSE",x=" ", title = "Scenario 2")+
  theme(
    axis.text.x.bottom = element_blank(),
    axis.text.y.left = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.ticks.x = element_blank(),
    plot.title = element_text( size=14,hjust = 0.5),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")
  )

plot.list.seg[[2]] <- RMSE_box
## plot the JSD
library(philentropy)
source("./simulation/plot/get_jsd.R")
nspot <- nrow(Results.dec.seg[[1]][[1]][[1]])
jsd = matrix(0,nrow = nspot, ncol = length(method_names))
colnames(jsd) = method_names
for (k in 1:length(method_names)) {

  jsd_mat = matrix(0,nrow = nspot, ncol = niter)
  for (i in 1:niter) {
    temp1 <- Results.dec.seg[[i]][[1]][[k]][,cell_type_seg]
    temp2 <- sim_pbac_seg_muraro_wang$spot_data[[i]]$sim_spots_pro[,cell_type_seg]
    jsd_mat[,i] = get_jsd(temp1, temp2)
  }
  jsd [,k] = apply(jsd_mat,1,mean)
}

quanjsd = apply(na.omit(jsd),2, quantile)
rownames(quanjsd)
colnames( quanjsd) = method_names
sort(quanjsd[3,] )
meanjsd  = apply(jsd,2, mean)
names(meanjsd) = method_names
sort(meanjsd)

JSD_en <- data.frame(JSD = c(jsd),
                     lable = factor(rep(method_names, each = nspot),
                                    levels = method_names))

jsd.ggplot <-ggplot(JSD_en, aes(x =lable, y = JSD))+
  geom_violin(col="black") + aes(fill = lable)+
  labs(title="Scenario 2", x="Method",
       y="JSD") + theme_classic()+
  theme(
    axis.text.x.bottom = element_text(size = 14,hjust = 1,angle =45),
    axis.text.y.left = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.ticks.x = element_blank(),
    plot.title = element_text( size=14,hjust = 0.5),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")
  )
# jsd.ggplot
plot.list.seg[[3]] <- jsd.ggplot
names(plot.list.seg) <- c("seg_PCC", "seg_RMSE", "seg_JSD")
# save(plot.list.seg, file = "ggplot.seg.evaluation.RData")

# patchwork::wrap_plots(error_box, RMSE_box, jsd.ggplot, nrow = 3, ncol = 1)
pfig2_temp <- c(plot.list.seg[1], plot.list.seg[2], plot.list.seg[3])
pfig2 <- cowplot::plot_grid(  plotlist = pfig2_temp, labels = letters[1:3], nrow = 3, label_size = 20)
print(pfig2)
# ggsave(pfig2, filename = "evaluation.seg.pdf", width = 12, height = 14)
# ggsave(pfig2, filename = "evaluation.seg.eps", width = 12, height = 14)

###### plot dot
w.mean <- Reduce("+", list_w_ensemble_seg)/length(list_w_ensemble_seg)
correlation.mean <- apply(R_square_seg,2, mean)[1:length(w.mean)]

data_bar = data.frame(Method = factor(names(w.mean),levels = names(w.mean)),
                      weight_est = w.mean,
                      PCC = correlation.mean)
library(ggpubr)
library(ggrepel)
p2 = ggplot(data_bar, aes(x = PCC , y = weight_est, color = factor(Method))) +
  geom_point(size = 3) +
  theme_classic() +
  labs(y="Weights assigned to \n individual methods",
       x="PCC of individual methods",title = "Segerstolpe") +
  theme(axis.text.x = element_text(size = 12, hjust = 1),
        axis.text.y.left = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        panel.border = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 14)) +
  geom_text_repel(aes(label = Method), size = 3.5) +
  stat_cor(method='pearson', cor.coef.name = "tau", size = 4,
           label.x = 0.2, label.y = 0.075,color='red') +
  stat_cor(method='spearman', cor.coef.name = 'rho', size = 4,
           label.x = 0.2, label.y = .065,color='red')
print(p2)
# ggsave(p2, filename = "wplot.pdf", width = 6, height = 5)
# ggsave(p2, filename = "wplot.eps", width = 6, height = 5)
#
# plot.list.w.seg <- p2
# save(plot.list.w.seg, file = "ggplot.seg.weight.RData")


