# Script to calculate and plot posterior occupancy distributions for 
# Old Growth and Secondary Forest over five years. It integrates MCMC 
# outputs from BirdNET and PROTAX Sound to estimate species richness 
# with 95% credible intervals. 

# Read data; load functions and packages
compdata <- readRDS("Utils/compdata.rds")
source("Functions/GetPsiMat.R")
source("Functions/GetPsisArray.R")
library(ggplot2)

# --------
# BirdNET
# --------
outbn <- readRDS("Outputs/ManausMSOD_MCMCchains.rds")

# Compute basic parameters from JAGS output
n.iter <- length(outbn$sims.list$beta.t)
n.year <- length(outbn$mean$mu.alpha.lpp.year)
n.spec <- length(outbn$mean$alpha.lpsi1)

# Old growth
psi1.og <- get_psi1_mat(outbn, "OG")
psis.og <- get_psis_arr(outbn, "OG")

# Create an array to bind all psis posterior information
allpsi.og <- array(NA, dim=c(n.spec,5,n.iter))
allpsi.og[,1,] <- psi1.og
allpsi.og[,2:5,] <- psis.og

# posterior summary
sum_psi_og <- apply(allpsi.og, c(2,3), sum)
mean_og <- apply(sum_psi_og, 1, mean)
lb_og <- apply(sum_psi_og, 1, quantile, 0.025)
ub_og <- apply(sum_psi_og, 1, quantile, 0.975)

# Secondary forest
psi1.sf <- get_psi1_mat(outbn,"SF")
psis.sf <- get_psis_arr(outbn,"SF")

allpsi.sf <- array(NA, dim=c(n.spec,5,n.iter))
allpsi.sf[,1,] <- psi1.sf
allpsi.sf[,2:5,] <- psis.sf

sum_psi_sf <- apply(allpsi.sf, c(2,3), sum)
mean_sf <- apply(sum_psi_sf, 1, mean)
lb_sf <- apply(sum_psi_sf, 1, quantile, 0.025)
ub_sf <- apply(sum_psi_sf, 1, quantile, 0.975)

# Data frame containing psis information for each forest type
psis_df <- data.frame(
  year = rep(1:5, 2),
  mean = c(mean_og, mean_sf),
  lb = c(lb_og, lb_sf),
  ub = c(ub_og, ub_sf),
  forest = rep(c("OG", "SF"), each = 5),
  dataset = "BirdNET"
)

# Plot BirdNET output
ggplot(psis_df) +
  geom_errorbar(aes(x = year, ymin = lb, ymax = ub, group = forest), 
                width=0.35, linewidth=0.3,
                position=position_dodge(0.4)) +
  geom_point(aes(x=year,y=mean,fill=forest), 
             size=3, shape=21,
             position=position_dodge(0.4)) +
  ylim(c(20, 60)) +
  #geom_line(aes(x = year, y = mean, group = forest, color = forest),
  #          position = position_dodge(0.4),
  #          linewidth = 0.4) +
  scale_fill_manual(name="Forest\ntype", values = c("black", "white")) +
  #scale_color_manual(values = c("OG" = "black",
  #                              "SF" = "grey60")) +
  #scale_y_log10(limits = c(20, 60)) + # y in log scale & y > 1
  labs(x = "Year",
       y = "Expected species richness per site",
       title = "Species richness per forest type and year - BirdNET data") +
  theme_bw() +
    theme(aspect.ratio = 1,
          axis.text.x = element_text(size=10,vjust=-1),
          axis.title.x = element_text(size = 12, vjust = -2),
          axis.text.y = element_text(size=9),
          axis.title.y = element_text(size=12,vjust=5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank())  
## clean up
rm(outbn)


# --------
# Montpellier
# --------
# Load Montpellier output: file named out16X
# load("C:/Users/gferr/OneDrive/Projects/NSpsOccDynManaus/MontpellierPaper/Outputs/output_4_WriteRunJAGSModel.RData")
load("OriginalData/output_4_WriteRunJAGSModel.RData")


# Compute parameters
n.iter.m <- length(out16X$sims.list$beta.t)
n.year.m <- length(out16X$mean$mu.alpha.lpp.year)
n.spec.m <- length(out16X$mean$alpha.lpsi1)

# Old growth
psi1.og.m <- get_psi1_mat(out16X, "OG")
psis.og.m <- get_psis_arr(out16X, "OG")

allpsi.og.m <- array(NA, dim=c(n.spec.m,5,n.iter.m))
allpsi.og.m[,1,] <- psi1.og.m
allpsi.og.m[,2:5,] <- psis.og.m

sum_psi_og.m <- apply(allpsi.og.m, c(2,3), sum)
mean_og.m <- apply(sum_psi_og.m, 1, mean)
lb_og.m <- apply(sum_psi_og.m, 1, quantile, 0.025)
ub_og.m <- apply(sum_psi_og.m, 1, quantile, 0.975)

# Secondary forest
psi1.sf.m <- get_psi1_mat(out16X, "SF")
psis.sf.m <- get_psis_arr(out16X, "SF")

allpsi.sf.m <- array(NA, dim=c(n.spec.m,5,n.iter.m))
allpsi.sf.m[,1,] <- psi1.sf.m
allpsi.sf.m[,2:5,] <- psis.sf.m

sum_psi_sf.m <- apply(allpsi.sf.m, c(2,3), sum)
mean_sf.m <- apply(sum_psi_sf.m, 1, mean)
lb_sf.m <- apply(sum_psi_sf.m, 1, quantile, 0.025)
ub_sf.m <- apply(sum_psi_sf.m, 1, quantile, 0.975)

psis_df_m <- data.frame(
  year = rep(1:5, 2),
  mean = c(mean_og.m, mean_sf.m),
  lb = c(lb_og.m, lb_sf.m),
  ub = c(ub_og.m, ub_sf.m),
  forest = rep(c("OG", "SF"), each = 5),
  dataset = "PROTAX Sound"
)

ggplot(psis_df_m) +
  geom_errorbar(aes(x = year, ymin = lb, ymax = ub, group = forest), 
                width=0.35, linewidth=0.3,
                position=position_dodge(0.4)) +
  geom_point(aes(x=year,y=mean,fill=forest), 
             size=3, shape=21,
             position=position_dodge(0.4)) +
  ylim(c(10, 60)) +
  scale_fill_manual(name="Forest\ntype", values = c("black", "white")) +
  labs(x = "Year",
       y = "Expected species richness per site",
       title = "Species richness per forest type and year - Montpellier data") +
  theme_bw() +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size=10,vjust=-1),
        axis.title.x = element_text(size = 12, vjust = -2),
        axis.text.y = element_text(size=9),
        axis.title.y = element_text(size=12,vjust=5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())  

# --------
# Bind the two projects data
# --------
psis_all <- rbind(psis_df, psis_df_m)

# Plot them together
ggplot(psis_all) +
  geom_errorbar(aes(x = year, ymin = lb, ymax = ub, group = forest),
                width = 0.35, linewidth = 0.3,
                position = position_dodge(0.4)) +
  geom_point(aes(x = year, y = mean, fill = forest),
             size = 3, shape = 21,
             position = position_dodge(0.4)) +
  ylim(c(10, 60)) +
  scale_fill_manual(name="Forest\ntype", values = c("black", "white")) +
  labs(x = "Year",
       y = "Expected species richness per site",
       title = "Species richness - BirdNET & PROTAX Sound") +
  facet_wrap(~dataset) + # facet_wrap: split by dataset
  theme_bw() +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size=10, vjust=-1),
        axis.title.x = element_text(size = 12, vjust = -2),
        axis.text.y = element_text(size=9),
        axis.title.y = element_text(size=12, vjust=5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  theme(plot.margin = margin(5.5, 5.5, 5.5, 20))

# Save it
ggsave("Outputs/species_richness_plot.png",
       width = 9,
       height = 6,
       dpi = 300,
       units = "in",
       limitsize = FALSE)
