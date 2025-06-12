library(causalDisco)
library(xtable)
library(ggplot2)

source("R/NCtools.R")
source("R/misc.R")

#################################################################################
## Computations for Section 2 (5 node DAG example) ##############################
#################################################################################

d <- 5

# "True" 5-node DAG
set.seed(1)
rantrue <- ncDAG(d, nedges = 8, permute = TRUE)
pcalg::plot(as.graphNEL(rantrue))

# "Estimated" 5-node DAG
set.seed(2) 
ranest <- ncDAG(d, nedges = 7, permute = TRUE)
pcalg::plot(as.graphNEL(ranest))

# Table 1 (adj. confusion matrix)
table1 <- confusion(ranest, rantrue)
precision(table1)
recall(table1)

# Simulation-based results
nrep <- 1000
precs_adj <- rep(NA, nrep)
recs_adj <- rep(NA, nrep)

set.seed(1)
for (i in 1:nrep) {
  ranest <- ncDAG(d, nedges = 7)
  thisconf_adj <- confusion(ranest, rantrue)
  precs_adj[i] <- precision(thisconf_adj)  
  recs_adj[i] <- recall(thisconf_adj)
}

median(precs_adj)
median(recs_adj)

#################################################################################
## Computations for Section 4 (Adjacency precision and recall) ##################
#################################################################################

# Example: Expected precision and recall for a dense 5 node DAG skeleton ########

nc_adj_precision(mest = 7, mtrue = 8, d = 5)
nc_adj_recall(mest = 7, mtrue = 8, d = 5)

# Example: F1 scores for a 5 node DAG with varying density ######################
d = 5
m_max <- maxnedges(d)
allnedges <- 1:m_max

d5f1res <- data.frame(m_true = rep(allnedges, m_max),
                      m_est = rep(allnedges, each = m_max),
                      d = 5)

d5f1res$adj_f1 <- nc_adj_f1(d5f1res$m_est, d5f1res$m_true, 5)$expectation

ggplot(d5f1res, aes(x = m_est, y = adj_f1, group = m_true, 
                    col = factor(m_true))) +  
  geom_point(size = 2) +
  geom_line() +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(breaks = seq(0,1,0.1)) + 
  xlab(bquote("Number of estimated edges"~(m[est]))) + 
  ylab("Expected adjacency F1") +
  guides(color = guide_legend(position = "bottom", nrow = 1)) + 
  scale_color_manual(
    name = expression(m[true]),
    values = c("#4F76A3", "#4C9A78", "#B0B836", "#D17C7F", "#B56E50",  
               "#9574A2", "#C1A193", "#6F6F6F", "#E0AE5C", "#A1426C")) + 
  theme_light()

ggsave(filename = "figures/d5anres-f1.pdf", height = 4, width = 6)

#################################################################################
## Computations for Section 5 (Metropolit application) ##########################
#################################################################################

# Metropolit application confuson matrix
metro_conf <- list(tp = 10, fp = 20, fn = 20, tn = 181)

# Compute p-value for overall test of skeleton fit
skelfit.test(metro_conf)

#################################################################################
## Computations for Section 6.1 (PC evaluation) #################################
#################################################################################

allres_pc_mtrue30 <- data.frame(mtrue = 30,
                            nedges = rep(NA, 1000),
                            shd = rep(NA, 1000),
                            vstructp = rep(NA, 1000),
                            ori_prec = rep(NA, 1000),
                            ori_rec = rep(NA, 1000),
                            adj_prec = rep(NA, 1000),
                            adj_rec = rep(NA, 1000),
                            sid_lower = rep(NA, 1000),
                            sid_upper = rep(NA, 1000))

allres_nc_mtrue30 <- allres_pc_mtrue30
allres_pc_mtrue30$method <- "pc"
allres_nc_mtrue30$method <- "nc"

allres_pc_mtrue15 <- allres_pc_mtrue30
allres_pc_mtrue15$mtrue <- 15
allres_nc_mtrue15 <- allres_nc_mtrue30
allres_nc_mtrue15$mtrue <- 15

# mtrue = 30 case ###############################################################
set.seed(1234)
rantrues_mtrue30 <- list()

for (i in 1:1000) {
  d <- 10
  m_max <- maxnedges(d)
  
  rantrue <- ncDAG(d, nedges = 30, permute = FALSE)
  rantrues_mtrue30 <- c(rantrues_mtrue30, list(rantrue))
  simdata <- simGausFromDAG(rantrue, n = 1000)
  
  pc_res <- amat(pc(simdata, sparsity = 0.05, test = corTest))
  
  sid <- SID::structIntervDist(rantrue, pc_res)
  confus_dir <- confusion(pc_res, rantrue, type = "dir")
  confus_adj <- confusion(pc_res, rantrue, type = "adj")
  
  allres_pc_mtrue30[i, "nedges"] <- nedges(pc_res)
  allres_pc_mtrue30[i, "shd"] <- shd(pc_res, rantrue)
  allres_pc_mtrue30[i, "sid_lower"] <- sid$sidLowerBound
  allres_pc_mtrue30[i, "sid_upper"] <- sid$sidUpperBound
  allres_pc_mtrue30[i, "ori_prec"] <- precision(confus_dir)
  allres_pc_mtrue30[i, "ori_rec"] <- recall(confus_dir)
  allres_pc_mtrue30[i, "adj_prec"] <- precision(confus_adj)
  allres_pc_mtrue30[i, "adj_rec"] <- recall(confus_adj)
  allres_pc_mtrue30[i, "vstructp"] <- percentVstruct(pc_res, rantrue)
  
}

set.seed(2345)

for (i in 1:1000) {
  rantrue <- rantrues_mtrue30[[i]]
  ranest <- ncCPDAG(d, nedges = sample(allres_pc_mtrue30$nedges, size = 1))
  
  sid <- SID::structIntervDist(rantrue, ranest)
  confus_dir <- confusion(ranest, rantrue, type = "dir")
  confus_adj <- confusion(ranest, rantrue, type = "adj")
  
  allres_nc_mtrue30[i, "nedges"] <- nedges(ranest)
  allres_nc_mtrue30[i, "shd"] <- shd(ranest, rantrue)
  allres_nc_mtrue30[i, "sid_lower"] <- sid$sidLowerBound
  allres_nc_mtrue30[i, "sid_upper"] <- sid$sidUpperBound
  allres_nc_mtrue30[i, "ori_prec"] <- precision(confus_dir)
  allres_nc_mtrue30[i, "ori_rec"] <- recall(confus_dir)
  allres_nc_mtrue30[i, "adj_prec"] <- precision(confus_adj)
  allres_nc_mtrue30[i, "adj_rec"] <- recall(confus_adj)
  allres_nc_mtrue30[i, "vstructp"] <- percentVstruct(ranest, rantrue)
  
}

# Range and mean of nedges
range(allres_nc_mtrue30$nedges)
mean(allres_nc_mtrue30$nedges)

# Organize results and export for Latex table
res_6.1_mtrue30 <- data.frame(
  metric = c("SHD", 
             "Adjacency precision",
             "Adjacency recall",
             "Orientation precision",
             "Orientation recall",
             "Perc. recovered v-structures",
             "SID (lower bound)",
             "SID (upper bound)"),
  PC_mean = format(c(mean(allres_pc_mtrue30$shd),
                     mean(allres_pc_mtrue30$adj_prec),
                     mean(allres_pc_mtrue30$adj_rec),
                     mean(allres_pc_mtrue30$ori_prec),
                     mean(allres_pc_mtrue30$ori_rec),
                     mean(allres_pc_mtrue30$vstructp),
                     mean(allres_pc_mtrue30$sid_lower),
                     mean(allres_pc_mtrue30$sid_upper)),
                   digits = 1),
  PC_CI = c(ci_latex_info(allres_pc_mtrue30$shd),
            ci_latex_info(allres_pc_mtrue30$adj_prec),
            ci_latex_info(allres_pc_mtrue30$adj_rec),
            ci_latex_info(allres_pc_mtrue30$ori_prec),
            ci_latex_info(allres_pc_mtrue30$ori_rec),
            ci_latex_info(allres_pc_mtrue30$vstructp),
            ci_latex_info(allres_pc_mtrue30$sid_lower),
            ci_latex_info(allres_pc_mtrue30$sid_upper)),
  NC_mean = format(c(mean(allres_nc_mtrue30$shd),
                     mean(allres_nc_mtrue30$adj_prec),
                     mean(allres_nc_mtrue30$adj_rec),
                     mean(allres_nc_mtrue30$ori_prec),
                     mean(allres_nc_mtrue30$ori_rec),
                     mean(allres_nc_mtrue30$vstructp),
                     mean(allres_nc_mtrue30$sid_lower),
                     mean(allres_nc_mtrue30$sid_upper)),
                   digits = 1),
  NC_CI = c(ci_latex_info(allres_nc_mtrue30$shd),
            ci_latex_info(allres_nc_mtrue30$adj_prec),
            ci_latex_info(allres_nc_mtrue30$adj_rec),
            ci_latex_info(allres_nc_mtrue30$ori_prec),
            ci_latex_info(allres_nc_mtrue30$ori_rec),
            ci_latex_info(allres_nc_mtrue30$vstructp),
            ci_latex_info(allres_nc_mtrue30$sid_lower),
            ci_latex_info(allres_nc_mtrue30$sid_upper)),
  p = format(c(mean(allres_nc_mtrue30$shd <= allres_pc_mtrue30$shd),
               mean(allres_nc_mtrue30$adj_prec >= allres_pc_mtrue30$adj_prec),
               mean(allres_nc_mtrue30$adj_rec >= allres_pc_mtrue30$adj_rec),
               mean(allres_nc_mtrue30$ori_prec >= allres_pc_mtrue30$ori_prec),
               mean(allres_nc_mtrue30$ori_rec >= allres_pc_mtrue30$ori_rec),
               mean(allres_nc_mtrue30$vstructp >= allres_pc_mtrue30$vstructp),
               mean(allres_nc_mtrue30$sid_lower <= allres_pc_mtrue30$sid_lower),
               mean(allres_nc_mtrue30$sid_upper <= allres_pc_mtrue30$sid_upper)),
             digits = 3)
  )

print(xtable(res_6.1_mtrue30), include.rownames = FALSE,
      sanitize.text.function = function(x) {x})


# mtrue = 15 case ###############################################################
set.seed(12345)
rantrues_mtrue15 <- list()

for (i in 1:1000) {
  d <- 10
  m_max <- maxnedges(d)
  
  rantrue <- ncDAG(d, nedges = 15, permute = FALSE)
  rantrues_mtrue15 <- c(rantrues_mtrue15, list(rantrue))
  simdata <- simGausFromDAG(rantrue, n = 1000)
  
  pc_res <- amat(pc(simdata, sparsity = 0.05, test = corTest))
  
  sid <- SID::structIntervDist(rantrue, pc_res)
  confus_dir <- confusion(pc_res, rantrue, type = "dir")
  confus_adj <- confusion(pc_res, rantrue, type = "adj")
  
  allres_pc_mtrue15[i, "nedges"] <- nedges(pc_res)
  allres_pc_mtrue15[i, "shd"] <- shd(pc_res, rantrue)
  allres_pc_mtrue15[i, "sid_lower"] <- sid$sidLowerBound
  allres_pc_mtrue15[i, "sid_upper"] <- sid$sidUpperBound
  allres_pc_mtrue15[i, "ori_prec"] <- precision(confus_dir)
  allres_pc_mtrue15[i, "ori_rec"] <- recall(confus_dir)
  allres_pc_mtrue15[i, "adj_prec"] <- precision(confus_adj)
  allres_pc_mtrue15[i, "adj_rec"] <- recall(confus_adj)
  allres_pc_mtrue15[i, "vstructp"] <- percentVstruct(pc_res, rantrue)
  
}

set.seed(23456)

for (i in 1:1000) {
  rantrue <- rantrues_mtrue15[[i]]
  ranest <- ncCPDAG(d, nedges = sample(allres_pc_mtrue15$nedges, size = 1))
  
  sid <- SID::structIntervDist(rantrue, ranest)
  confus_dir <- confusion(ranest, rantrue, type = "dir")
  confus_adj <- confusion(ranest, rantrue, type = "adj")
  
  allres_nc_mtrue15[i, "nedges"] <- nedges(ranest)
  allres_nc_mtrue15[i, "shd"] <- shd(ranest, rantrue)
  allres_nc_mtrue15[i, "sid_lower"] <- sid$sidLowerBound
  allres_nc_mtrue15[i, "sid_upper"] <- sid$sidUpperBound
  allres_nc_mtrue15[i, "ori_prec"] <- precision(confus_dir)
  allres_nc_mtrue15[i, "ori_rec"] <- recall(confus_dir)
  allres_nc_mtrue15[i, "adj_prec"] <- precision(confus_adj)
  allres_nc_mtrue15[i, "adj_rec"] <- recall(confus_adj)
  allres_nc_mtrue15[i, "vstructp"] <- percentVstruct(ranest, rantrue)
  
}

# Range and mean of nedges
range(allres_nc_mtrue15$nedges)
mean(allres_nc_mtrue15$nedges)

# Organize results and export for Latex table
res_6.1_mtrue15 <- data.frame(
  metric = c("SHD", 
             "Adjacency precision",
             "Adjacency recall",
             "Orientation precision",
             "Orientation recall",
             "Perc. recovered v-structures",
             "SID (lower bound)","SID (upper bound)"),
  PC_mean = format(c(mean(allres_pc_mtrue15$shd),
                     mean(allres_pc_mtrue15$adj_prec),
                     mean(allres_pc_mtrue15$adj_rec),
                     mean(allres_pc_mtrue15$ori_prec),
                     mean(allres_pc_mtrue15$ori_rec),
                     mean(allres_pc_mtrue15$vstructp),
                     mean(allres_pc_mtrue15$sid_lower),
                     mean(allres_pc_mtrue15$sid_upper)),
                   digits = 1),
  PC_CI = c(ci_latex_info(allres_pc_mtrue15$shd),
            ci_latex_info(allres_pc_mtrue15$adj_prec),
            ci_latex_info(allres_pc_mtrue15$adj_rec),
            ci_latex_info(allres_pc_mtrue15$ori_prec),
            ci_latex_info(allres_pc_mtrue15$ori_rec),
            ci_latex_info(allres_pc_mtrue15$vstructp),
            ci_latex_info(allres_pc_mtrue15$sid_lower),
            ci_latex_info(allres_pc_mtrue15$sid_upper)),
  NC_mean = format(c(mean(allres_nc_mtrue15$shd),
                     mean(allres_nc_mtrue15$adj_prec),
                     mean(allres_nc_mtrue15$adj_rec),
                     mean(allres_nc_mtrue15$ori_prec),
                     mean(allres_nc_mtrue15$ori_rec),
                     mean(allres_nc_mtrue15$vstructp),
                     mean(allres_nc_mtrue15$sid_lower),
                     mean(allres_nc_mtrue15$sid_upper)),
                   digits = 1),
  NC_CI = c(ci_latex_info(allres_nc_mtrue15$shd),
            ci_latex_info(allres_nc_mtrue15$adj_prec),
            ci_latex_info(allres_nc_mtrue15$adj_rec),
            ci_latex_info(allres_nc_mtrue15$ori_prec),
            ci_latex_info(allres_nc_mtrue15$ori_rec),
            ci_latex_info(allres_nc_mtrue15$vstructp),
            ci_latex_info(allres_nc_mtrue15$sid_lower),
            ci_latex_info(allres_nc_mtrue15$sid_upper)),
  p = format(c(mean(allres_nc_mtrue15$shd <= allres_pc_mtrue15$shd),
               mean(allres_nc_mtrue15$adj_prec >= allres_pc_mtrue15$adj_prec),
               mean(allres_nc_mtrue15$adj_rec >= allres_pc_mtrue15$adj_rec),
               mean(allres_nc_mtrue15$ori_prec >= allres_pc_mtrue15$ori_prec),
               mean(allres_nc_mtrue15$ori_rec >= allres_pc_mtrue15$ori_rec),
               mean(allres_nc_mtrue15$vstructp >= allres_pc_mtrue15$vstructp),
               mean(allres_nc_mtrue15$sid_lower <= allres_pc_mtrue15$sid_lower),
               mean(allres_nc_mtrue15$sid_upper <= allres_pc_mtrue15$sid_upper)),
             digits = 3)
)

print(xtable(res_6.1_mtrue15), include.rownames = FALSE,
      sanitize.text.function = function(x) {x})



#################################################################################
## Computations for Section 6.2 (Sachs application) #############################
#################################################################################

# Import ground truth graph
sachs_truth <- as.matrix(read.table("data/sachs.2005.truth.cpdagadj.txt", 
                                    header = TRUE))

# Import Tetrad results for BOSS, LiNGAM, GES and PC
# Tetrad session can be found in "tetrad/sachs-application.tet"
boss_res <- as.matrix(read.table("data/BOSS.amat.cpag-sachs.txt", 
                                 header = TRUE))
lingam_res <- as.matrix(read.table("data/lingam.amat.cpag-sachs.txt", 
                                   header = TRUE))
ges_res <- as.matrix(read.table("data/ges.amat.cpag-sachs.txt", 
                                header = TRUE))
pc_res <- as.matrix(read.table("data/pc.amat.cpag-sachs.txt", 
                               header = TRUE))

# Compute SHDs and count estimated number of edges
shd_obs_boss <- shd(boss_res, sachs_truth)
shd_obs_lingam <- shd(lingam_res, sachs_truth)
shd_obs_ges <- shd(ges_res, sachs_truth)
shd_obs_pc <- shd(pc_res, sachs_truth)

nedges_boss <- nedges(boss_res)
nedges_lingam <- nedges(lingam_res)
nedges_ges <- nedges(ges_res)
nedges_pc <- nedges(pc_res)

shd_obs_notears <- 22 #from notears paper
nedges_notears <- 16#from notears paper

# Compute neg. control SHDs for each algorithm
nrep <- 1000
shds_boss <- rep(NA, nrep)
shds_ges <- rep(NA, nrep)
shds_lingam <- rep(NA, nrep)
shds_pc <- rep(NA, nrep)
shds_notears <- rep(NA, nrep)

set.seed(1)
for (i in 1:nrep) {
  shds_notears[i] <- shd(ncDAG(11, nedges = nedges_notears), sachs_truth)
  shds_boss[i] <- shd(ncCPDAG(11, nedges = nedges_boss), sachs_truth)
  shds_ges[i] <- shd(ncCPDAG(11, nedges = nedges_ges), sachs_truth)
  shds_lingam[i] <- shd(ncDAG(11, nedges = nedges_lingam), sachs_truth)
  shds_pc[i] <- shd(ncCPDAG(11, nedges = nedges_pc), sachs_truth)
}

# Combine results in table
res_6.2 <- data.frame(algo = c("NOTEARS", "BOSS", "LiNGAM", "GES", "PC"),
                      obs_SHD = c(shd_obs_notears, shd_obs_boss, 
                                  shd_obs_lingam, shd_obs_ges,
                                  shd_obs_pc),
                      obs_mest = c(nedges_notears, nedges_boss,
                                   nedges_lingam, nedges_ges, nedges_pc),
                      nc_SHD = c(mean(shds_notears), mean(shds_boss),
                                  mean(shds_lingam), mean(shds_ges), 
                                  mean(shds_pc)),
                      p = format(c(mean(shds_notears <= shd_obs_notears),
                            mean(shds_boss <= shd_obs_boss),
                            mean(shds_lingam <= shd_obs_lingam),
                            mean(shds_ges <= shd_obs_ges),
                            mean(shds_pc <= shd_obs_pc)), digits = 3)
                      )

# order by p-values
res_6.2 <- res_6.2[order(res_6.2$p),]

#export to Latex
print(xtable(res_6.2), include.rownames = FALSE)

      