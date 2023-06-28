source("functions_max_matching.R")
library(cowplot)
library(ggbeeswarm)
library(tidyverse)
library(RImageJROI)
library(igraph)

# read ROI file
x <- read.ijroi("TFC_ROI_63Pix-120um.roi")

## Redefine ROI coordinates manually to cover more precisely target area
x$xrange <- c(330, 500)
x$yrange <- c(0,450)
x$coords[,1] <- c(320,320,450,450-320)
x$coords[,2] <- c(0,0,450,450)

c1_list_of_files <- list.files("tracking_data/gfp-tracks", pattern="*.csv")
# create a plots folder, if it doesnt exist within the project
dir.create(file.path(getwd(), "/plots"), showWarnings = FALSE)

# generate plots of scatterplots with ROI overlaid
for (c1_file in c1_list_of_files) {
  # find corresponding file in mCherry folder
  c3_file <- gsub("^C1","C3", c1_file )
  
  f1 <- paste0("tracking_data/gfp-tracks/", c1_file)
  f2 <- paste0("tracking_data/mcherry-tracks/", c3_file)
  
  # Load both track files
  d1 <- read.csv( f1 )
  d2 <- read.csv( f2 )
  
  name <-paste0("plots/scatter_C3_",substring(c1_file,1,nchar(c1_file)-4),".pdf")
  pdf(name, width=5, height=5.5, pointsize = 9, useDingbats=FALSE)
  plot(d1$POSITION_X,d1$POSITION_Y, xlab = "POSITION_X", ylab="POSITION_Y", 
       main = paste0("ROI over points-", "C3_",substring(c1_file,1,nchar(c1_file)-4) ),
       cex.main=0.7,
       cex.lab=0.8)
  points(d2$POSITION_X,d2$POSITION_Y, pch=19, cex=.5)
  rect( x$xrange[1], x$yrange[1], x$xrange[2], x$yrange[2], border=2, lwd=3)
  #plot(x, border = "red", add = TRUE)
  dev.off()
  
}

mix_02h <- c()
nt_02h <- c()
rescue_02h <- c()

mix_24h <- c() 
nt_24h <- c()
rescue_24h <- c()

mix_46h <- c() 
nt_46h <- c()
rescue_46h <- c()

mix_68h <- c() 
nt_68h <- c()
rescue_68h <- c()


# filter data and only focus on the ROI
for (c1_file in c1_list_of_files) {
  # find corresponding file in mCherry folder
  c3_file <- gsub("^C1","C3", c1_file )
  
  f1 <- paste0("tracking_data/gfp-tracks/", c1_file)
  f2 <- paste0("tracking_data/mcherry-tracks/", c3_file)
  
  # Load both track files
  d1 <- read.csv( f1 )
  d2 <- read.csv( f2 )
  
  # filter according to ROI coordinates
  d1 <- d1[d1$POSITION_X>= x$xrange[1],]
  d1 <- d1[d1$POSITION_X<=x$xrange[2],]
  d1 <- d1[d1$POSITION_Y>=x$yrange[1],]
  d1 <- d1[d1$POSITION_Y<=x$yrange[2],]
  
  d2 <- d2[d2$POSITION_X>= x$xrange[1],]
  d2 <- d2[d2$POSITION_X<=x$xrange[2],]
  d2 <- d2[d2$POSITION_Y>=x$yrange[1],]
  d2 <- d2[d2$POSITION_Y<=x$yrange[2],]
  
  # note that min and max time are given in seconds
  print(c1_file)
  plt_02h <- extract_time_interval(d1, d2, min_time=0, max_time=7200, min_time_incl=TRUE, max_time_incl=TRUE,min_track_length=5, min_dist_cells=10)
  plt_24h <- extract_time_interval(d1, d2, min_time=7200, max_time=14400, min_time_incl=FALSE, max_time_incl=TRUE,min_track_length=5, min_dist_cells=10)
  plt_46h <- extract_time_interval(d1, d2, min_time=14400, max_time=21600, min_time_incl=FALSE, max_time_incl=TRUE,min_track_length=5, min_dist_cells=10)
  plt_68h <- extract_time_interval(d1, d2, min_time=21600, max_time=28800, min_time_incl=FALSE, max_time_incl=TRUE,min_track_length=5, min_dist_cells=10)
  print(plt_46h)
 
  if (grepl("MIX", c1_file, ignore.case=TRUE)) { 
    mix_02h <- append(mix_02h, plt_02h)
    mix_24h <- append(mix_24h, plt_24h)     
    mix_46h <- append(mix_46h, plt_46h)
    mix_68h <- append(mix_68h, plt_68h)
  } else if (grepl("NT", c1_file, ignore.case=TRUE)) {
    nt_02h <- append(nt_02h, plt_02h)
    nt_24h <- append(nt_24h, plt_24h)
    nt_46h <- append(nt_46h, plt_46h)
    nt_68h <- append(nt_68h, plt_68h)
  } else if (grepl("RESCUE", c1_file, ignore.case=TRUE)) {
    rescue_02h <- append(rescue_02h, plt_02h)
    rescue_24h <- append(rescue_24h, plt_24h)
    rescue_46h <- append(rescue_46h, plt_46h)
    rescue_68h <- append(rescue_68h, plt_68h)
  }
}


all_02h <- tibble(mean_int_len=append(append(rescue_02h,mix_02h),nt_02h),
                  type=c(rep("rescue",length(rescue_02h)), rep("mix",length(mix_02h)),rep("nt",length(nt_02h))),
                  hour_int=rep("0-2h", length(mean_int_len))
)

all_24h <- tibble(mean_int_len=append(append(rescue_24h,mix_24h),nt_24h),
                  type=c(rep("rescue",length(rescue_24h)), rep("mix",length(mix_24h)),rep("nt",length(nt_24h))),
                  hour_int=rep("2-4h", length(mean_int_len))
)


all_46h <- tibble(mean_int_len=append(append(rescue_46h,mix_46h),nt_46h),
                  type=c(rep("rescue",length(rescue_46h)), rep("mix",length(mix_46h)),rep("nt",length(nt_46h))),
                  hour_int=rep("4-6h", length(mean_int_len))
)


all_68h <- tibble(mean_int_len=append(append(rescue_68h,mix_68h),nt_68h),
                  type=c(rep("rescue",length(rescue_68h)), rep("mix",length(mix_68h)),rep("nt",length(nt_68h))),
                  hour_int=rep("6-8h", length(mean_int_len))
)

max_all_int <- ceiling(max(all_02h$mean_int_len,all_24h$mean_int_len,all_46h$mean_int_len,all_68h$mean_int_len))
# used for plots
break_value <- 2

plt_all_02h <- ggplot(all_02h, aes(x = factor(type, levels = c("nt", "mix", "rescue")), y = mean_int_len, colour=factor(type))) +
  geom_beeswarm()+
  theme_classic()+
  labs(x="",y="Interaction duration (min)", title="0-2h")+
  scale_y_continuous(limits = c(0,max_all_int), breaks = seq(0,max_all_int, by=break_value))+
  scale_x_discrete(labels = c('NT','PMM2 KD','Rescue')) +
  theme(legend.position="none")

plt_all_24h <- ggplot(all_24h, aes(x = factor(type, levels = c("nt", "mix", "rescue")), y = mean_int_len, colour=factor(type))) +
  geom_beeswarm()+
  theme_classic()+
  labs(x="",y="", title="2-4h")+
  scale_y_continuous(limits = c(0,max_all_int), breaks = seq(0,max_all_int, by=break_value))+
  scale_x_discrete(labels = c('NT','PMM2 KD','Rescue')) +
  theme(legend.position="none")

plt_all_46h <- ggplot(all_46h, aes(x = factor(type, levels = c("nt", "mix", "rescue")), y = mean_int_len, colour=factor(type))) +
  geom_beeswarm()+
  theme_classic()+
  labs(x="",y="", title="4-6h")+
  scale_y_continuous(limits = c(0,max_all_int), breaks = seq(0,max_all_int, by=break_value))+
  scale_x_discrete(labels = c('NT','PMM2 KD','Rescue')) +
  theme(legend.position="none")

plt_all_68h <- ggplot(all_68h, aes(x = factor(type, levels = c("nt", "mix", "rescue")), y = mean_int_len, colour=factor(type))) +
  geom_beeswarm()+
  theme_classic()+
  labs(x="",y="", title="6-8h")+
  scale_y_continuous(limits = c(0,max_all_int), breaks = seq(0,max_all_int, by=break_value))+
  scale_x_discrete(labels = c('NT','PMM2 KD','Rescue')) +
  #guides(color = guide_legend(title = "Type"))
  theme(legend.position="none")

p_ls = list(plt_all_02h, plt_all_24h,plt_all_46h, plt_all_68h)

name <-paste0("plots/mean_edge_weight_per_movie_per_type.pdf")
pdf(name, width=10, height=4, pointsize = 9, useDingbats=FALSE) #units="in",res=600)
print(plot_grid(plotlist = p_ls, align = 'h', ncol = 4))
dev.off()

# plot mean interaction vs time (one plot for each type + one combined)

#nt
all_nt <- tibble(mean_int_len=append(append(append(nt_02h,nt_24h),nt_46h),nt_68h),
                 time_int=c(rep("0-2h",length(nt_02h)), rep("2-4h",length(nt_24h)),rep("4-6h",length(nt_46h)), rep("6-8h",length(nt_68h))),
                 
)

# add mean and sd for mean_int_len per time_int
all_nt$mean_value <- ave(all_nt$mean_int_len, all_nt$time_int, FUN = mean)
all_nt$sd_value <- ave(all_nt$mean_int_len, all_nt$time_int, FUN = sd)
# add a numeric column which allows for plotting the smooth line (otherwise it doesnt work with the factor)
all_nt$time_num <- as.numeric(gsub("([0-9]+)-.*", "\\1", all_nt$time_int))



all_mix <- tibble(mean_int_len=append(append(append(mix_02h,mix_24h),mix_46h),mix_68h),
                  time_int=c(rep("0-2h",length(mix_02h)), rep("2-4h",length(mix_24h)),rep("4-6h",length(mix_46h)), rep("6-8h",length(mix_68h))),
                  
)

all_mix$mean_value <- ave(all_mix$mean_int_len, all_mix$time_int, FUN = mean)
all_mix$sd_value <- ave(all_mix$mean_int_len, all_mix$time_int, FUN = sd)
# add a numeric column which allows for plotting the smooth line (otherwise it doesnt work with the factor)
all_mix$time_num <- as.numeric(gsub("([0-9]+)-.*", "\\1", all_mix$time_int))


all_rescue <- tibble(mean_int_len=append(append(append(rescue_02h,rescue_24h),rescue_46h),rescue_68h),
                     time_int=c(rep("0-2h",length(rescue_02h)), rep("2-4h",length(rescue_24h)),rep("4-6h",length(rescue_46h)), rep("6-8h",length(rescue_68h))),
                     
)

all_rescue$mean_value <- ave(all_rescue$mean_int_len, all_rescue$time_int, FUN = mean)
all_rescue$sd_value <- ave(all_rescue$mean_int_len, all_rescue$time_int, FUN = sd)
# add a numeric column which allows for plotting the smooth line (otherwise it doesnt work with the factor)
all_rescue$time_num <- as.numeric(gsub("([0-9]+)-.*", "\\1", all_rescue$time_int))


all_conditions <- data.frame(bind_rows(all_nt, all_mix, all_rescue))
all_conditions$type <- c(rep("nt", length(all_nt$mean_int_len)), rep("mix", length(all_mix$mean_int_len)), rep("rescue", length(all_rescue$mean_int_len)))

max_all_int = max(all_nt$mean_int_len, all_mix$mean_int_len, all_rescue$mean_int_len)

# plots with beeswarms, smooth line and se band
#nt
plt_all_nt <- ggplot(all_nt, aes(x = time_num, y = mean_int_len)) +
  geom_beeswarm(size=1) +
  geom_smooth(method = "loess", se = FALSE, color = "#F8766D", linewidth=1) +
  geom_ribbon(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), 
              data = all_nt, alpha = 0.2, fill = "grey") +
  theme_classic()+
  labs(x="",y="Interaction duration (min)", title="NT")+
  scale_y_continuous(limits = c(0,max_all_int), breaks = seq(0,max_all_int, by=break_value))+
  scale_x_continuous(labels = c('0-2h', '2-4h', '4-6h', '6-8h'), breaks = c(0, 2, 4, 6))+
  theme(plot.title = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10))

#mix
plt_all_mix <- ggplot(all_mix, aes(x = time_num, y = mean_int_len)) +
  geom_beeswarm(size=1) +
  geom_smooth(method = "loess", se = FALSE, color = "#00BA38", linewidth=1) +
  geom_ribbon(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), 
              data = all_mix, alpha = 0.2, fill = "grey") +
  theme_classic()+
  labs(x="",y="", title="PMM2 KD")+
  scale_y_continuous(limits = c(0,max_all_int), breaks = seq(0,max_all_int, by=break_value))+
  scale_x_continuous(labels = c('0-2h', '2-4h', '4-6h', '6-8h'), breaks = c(0, 2, 4, 6))+
  theme(plot.title = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10))
#rescue
plt_all_rescue <- ggplot(all_rescue, aes(x = time_num, y = mean_int_len)) +
  geom_beeswarm(size=1) +
  geom_smooth(method = "loess", se = FALSE, color = "#619CFF", linewidth=1) +
  geom_ribbon(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), 
              data = all_rescue, alpha = 0.2, fill = "grey") +
  theme_classic()+
  labs(x="",y="", title="Rescue")+
  scale_y_continuous(limits = c(0,max_all_int), breaks = seq(0,max_all_int, by=break_value))+
  scale_x_continuous(labels = c('0-2h', '2-4h', '4-6h', '6-8h'), breaks = c(0, 2, 4, 6))+
  theme(plot.title = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10))

all_conditions$type <- factor(all_conditions$type, levels = c("nt", "mix", "rescue"))

plt_all <- ggplot(all_conditions, aes(x = time_num, y = mean_int_len, interaction(type, time_num), color = type)) +
  #geom_beeswarm(size=1) +
  geom_smooth(method = "loess", se = FALSE, linewidth=1) +
  geom_ribbon(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value, fill=type), colour = NA, 
              alpha = 0.4) +
  theme_classic() +
  labs(x = "", y = "", title = "All conditions") +
  scale_y_continuous(limits = c(0, max_all_int), breaks = seq(0, max_all_int, by=break_value)) +
  scale_x_continuous(labels = c('0-2h', '2-4h', '4-6h', '6-8h'), breaks = c(0, 2, 4, 6))+
  scale_color_discrete(labels=c('NT', 'PMM2 KD', 'Rescue'))+
  theme(plot.title = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        legend.title = element_text(size=8),
        legend.text=element_text(size=8))+
  guides(color=guide_legend("type"), fill = "none")+
  guides(color = guide_legend(title = "Type"))

print(plt_all)

p_all = list(plt_all_nt, plt_all_mix,plt_all_rescue, plt_all)

name <-paste0("plots/mean_edge_weight_maxmatching_vs_time_per_cond_smooth.pdf")
pdf(name, width=9, height=2, pointsize = 5, useDingbats=FALSE) #units="in",res=600)
print(plot_grid(plotlist = p_all, rel_widths=c(0.6,0.6,0.6,1.1),align = 'h', ncol = 4))
dev.off()