#script to assess parsed diversitool outputs

setwd("/home/rebee/projects/cyclophosphamide-mice/nimagen2/DiversiTools/")

library(scales)
library(ggplot2)
library(tidyverse)
library(gggenes)
library(patchwork)
library(ggrepel)
library(reshape2)
library(ggpubr)
library(rstatix)
library(nortest)
library(RColorBrewer)
library(ggbreak)
library(ggsignif)
library(cowplot)

####read in metadata

metadata_sp <- read.csv("../nextclade/metadata3.csv")
metadata_sp <- metadata_sp %>% filter(seqName!="Sample_NEGATIVE_NEGATIVE_1a")
metadata_sp <- metadata_sp %>% filter(seqName!="Sample_NEGATIVE_NEGATIVE_2a")
metadata_sp <- metadata_sp %>% filter(seqName!="Sample_POSITIVE_POSITIVE_1b")
metadata_sp <- metadata_sp %>% filter(seqName!="Sample_POSITIVE_POSITIVE_2b")


### Tv Ts analysis mice swabs ####
#entropy/TvTs inputs 


diversitools_readin_tvts <- function(filepath, file_pattern, file_list) {
  
  temp <- list.files(filepath, pattern=file_pattern) # creates list variable  with file names as they are in the directory
  file_list <- list() # creates empty list to add to in the for loop
  
  for (i in 1:length(temp)) { # read in all data files as DFs, add Tv and Ts ratio column to all 
    
    file_list[[i]] <- read.delim(temp[i]) 
    file_list[[i]]$Coverage <- as.numeric(file_list[[i]]$Coverage)
    file_list[[i]]$CntTs <- as.numeric(file_list[[i]]$CntTs)
    file_list[[i]]$CntTv <- as.numeric(file_list[[i]]$CntTv)
    file_list[[i]]$TsRatio <- file_list[[i]]$CntTs / file_list[[i]]$Coverage
    file_list[[i]]$TvRatio <- file_list[[i]]$CntTv / file_list[[i]]$Coverage
    file_list[[i]]$ins_cnt <- as.numeric(file_list[[i]]$ins_cnt)
    file_list[[i]]$del_cnt <- as.numeric(file_list[[i]]$del_cnt)
    file_list[[i]]$insRatio <- file_list[[i]]$ins_cnt / file_list[[i]]$Coverage
    file_list[[i]]$delRatio <- file_list[[i]]$del_cnt / file_list[[i]]$Coverage
    file_list[[i]][,6:16] <- sapply(file_list[[i]][,6:16], as.character)
    file_list[[i]]$Ncnt <- as.character(file_list[[i]]$Ncnt)
    
  }
  
  #these two lines will be specific to your divtool output file naming structure
  myfilenames <- gsub(".final_entropy.txt", "", temp) #removes common strings
  myfilenames <- gsub("_[^_]+$", "", temp) #removes the unique index strings after the first underscore, leaving just the kit_ids
  names(file_list) <- myfilenames
  return(file_list)
}
mut_list_sp <- diversitools_readin_tvts(filepath = "./", file_pattern = "entropy.txt", file_list = mut_list)

#changing divtools Sample column name and removing unnecessary strings from descriptor
for (i in 1:length(mut_list_sp)) {
  mut_list_sp[[i]]$Sample <- gsub(".final.bam", "", mut_list_sp[[i]]$Sample)
}


all_mut_sp <- bind_rows(mut_list_sp)
all_mut_sp <- merge(metadata_sp, all_mut_sp, by = "Sample")
all_mut_sp <- all_mut_sp %>% filter(Exclude!="Y")
all_mut_sp$Cohort <- factor(all_mut_sp$Cohort,levels=c("Vehicle", "Cyclophosphamide only", "Molnupiravir","Cyclophosphamide and Molnupiravir","Cyclophosphamide and Pfizer","Cyclophosphamide, Pfizer and Molnupiravir"))
melt_mut_sp <- reshape2::melt(all_mut_sp, na.rm=F, id.vars=c("Cohort", "Animal.label", "Sample.Type", "DPI", "Position", "Coverage"), measure.vars=c("TsRatio", "TvRatio"))# melt to long form
melt_mut_sp$variable <- as.factor(melt_mut_sp$variable)
melt_mut_sp$value <- as.numeric(melt_mut_sp$value)



#get Ts/Tv ratio per sample
meanTs_sp <- all_mut_sp %>% filter(Coverage >200) %>% group_by(Cohort, Animal.label, Sample.Type, DPI) %>% 
  summarize(Transition=mean(TsRatio, na.rm=T ), stdev=sd(TsRatio, na.rm=T)) 

meanTv_sp <- all_mut_sp %>% filter(Coverage >200) %>% group_by(Cohort, Animal.label, Sample.Type, DPI) %>% 
  summarize(Transversion=mean(TvRatio, na.rm=T ), stdev=sd(TvRatio, na.rm=T)) 

meanTsTv_sp <- merge(meanTs_sp, meanTv_sp, by =c("Cohort", "Animal.label", "Sample.Type", "DPI"))
meanTsTv_sp$TsTvRatio <- meanTsTv_sp$Transition / meanTsTv_sp$Transversion
melt_meanTsTv_sp <- reshape2::melt(meanTsTv_sp, na.rm=F, id.vars=c("Cohort", "Animal.label", "Sample.Type", "DPI"), measure.vars=c("TsTvRatio", "Transition", "Transversion"))


#### plot Ts/Tv ratios
library(ggsignif)
library(ggstatsplot)


# perform stats test for plot Cohorts faceted by DPI for throat swabs

stats1 <-melt_meanTsTv_sp %>%
  filter(variable=="TsTvRatio") %>% 
  filter(Sample.Type=="Swab") %>% 
  filter(DPI != "5") %>% 
  filter(DPI != "4") %>%
  rstatix::group_by(DPI) %>% 
  rstatix::wilcox_test(value~Cohort) %>% 
  adjust_pvalue(method="none") %>% 
  add_significance("p") %>%
  add_xy_position(x = "Cohort", dodge = 0.8) %>%
  remove_ns(col="p")

dodge <- position_dodge(width=1)

plot_tvtsRatio_stats_facet <- ggplot(melt_meanTsTv_sp %>% filter(variable=="TsTvRatio") %>% 
                                             filter(Sample.Type=="Swab") %>% filter(DPI != "5") %>% 
                                             filter(DPI != "4"), aes(x=Cohort, y=value, colour=Cohort))+ 
                          facet_grid(~DPI, scales="free") +
                          geom_boxplot(varwidth = TRUE, outlier.colour = "NA", alpha = 0.5)+
                          ylab("Transition/Transversion Ratio") +
                          xlab("Treatment") +
                          theme_bw() +
                          scale_colour_viridis_d()+
                          theme(axis.text.y = element_text(colour="black", size=12),
                                axis.title.y = element_text(size = 12, face="bold"),
                                axis.text.x = element_text(size=12, colour="black", angle=45, hjust=1),
                                axis.title.x = element_text(size = 12, face="bold"),
                                strip.background = element_rect(size=0.5, linetype="solid", fill="white"),
                                strip.text.x = element_text(size = 12, color = "black", face="bold"),
                                panel.grid.major = element_blank(), 
                                panel.grid.minor = element_blank(),
                                panel.background = element_blank(), 
                                axis.line = element_line(colour = "black"),
                                legend.position = "none") +
  stat_pvalue_manual(stats1, label = "p.signif")

        

plot_tvtsRatio_stats_facet


#write out 
#png for git 
ggsave(plot_tvtsRatio_stats_facet, filename="~/projects/cyclophosphamide-mice/nimagen2/figures2/CMP-tvtsratio_boxplot_facet-swabs-stats.png", device="png", dpi=300, width=8, height=8)
#tiff for pub
ggsave(plot_tvtsRatio_stats_facet, filename="~/projects/cyclophosphamide-mice/nimagen2/figures2/CMP-tvtsratio_boxplot_facet-swabs-stats.tiff", device="tiff", dpi=300, width=8, height=8)

# perform stats test for plot Cohorts faceted by tissue

stats2 <-melt_meanTsTv_sp %>%
  filter(variable=="TsTvRatio") %>% 
  filter(Sample.Type!="Swab") %>% 
  rstatix::group_by(Sample.Type) %>% 
  rstatix::wilcox_test(value~Cohort) %>% 
  adjust_pvalue(method="none") %>% 
  add_significance("p") %>%
  add_xy_position(x = "Cohort", dodge = 0.8) %>%
  remove_ns(col="p")



plot_tvtsRatio_stats_facet <- ggplot(melt_meanTsTv_sp %>% filter(variable=="TsTvRatio") %>% filter(Sample.Type!="Swab"), aes(x=Cohort, y=value, colour=Cohort))+ # filter out positions than have cov<10
  facet_grid(~Sample.Type, scales="free") +
  geom_boxplot(varwidth = TRUE, alpha = 0.5)+
  ylab("Transition/Transversion Ratio") +
  xlab("Treatment") +
  theme_bw() +
  scale_colour_viridis_d()+
  theme(axis.text.y = element_text(colour="black", size=12),
        axis.title.y = element_text(size = 12, face="bold"),
        axis.text.x = element_text(size=12, colour="black", angle=45, hjust=1),
        axis.title.x = element_text(size = 12, face="bold"),
        strip.background = element_rect(size=0.5, linetype="solid", fill="white"),
        strip.text.x = element_text(size = 12, color = "black", face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none") +
stat_pvalue_manual(stats2, label = "p.signif")


plot_tvtsRatio_stats_facet

#png for git 
ggsave(plot_tvtsRatio_stats_facetP_sp, filename="~/projects/cyclophosphamide-mice/nimagen2/figures2/CMP-tvtsratio_boxplot_facet-tissues-stats.png", device="png", dpi=300, width=8, height=8)
#tiff for pub
ggsave(plot_tvtsRatio_stats_facetP_sp, filename="~/projects/cyclophosphamide-mice/nimagen2/figures2/CMP-tvtsratio_boxplot_facet-tissues-stats.tiff", device="tiff", dpi=300, width=8, height=8)


#genome nucleotide position of proteins
nsps <- read.csv("../../bin/protein_position.csv") 
names(nsps)[1]<-'Position'
sapply(nsps,class) #checking class of positions
nsps$Position<-as.numeric(nsps$Position) # change pos to numeric

# sapply(all_mut, class)
all_mut<-merge(nsps,all_mut_sp)
all_mut$Acnt <- as.numeric(all_mut$Acnt)
all_mut$Ccnt <- as.numeric(all_mut$Ccnt)
all_mut$Gcnt <- as.numeric(all_mut$Gcnt)
all_mut$Tcnt <- as.numeric(all_mut$Tcnt)

#mutations#### makes indivdiual tibble for each nucleotide 
#G nucleotide 
all_mut_G <- subset(all_mut, RefBase == "G")
all_mut_G < as_tibble(all_mut_G)
#A nucleotide
all_mut_A <- subset(all_mut, RefBase == "A")
all_mut_A < as_tibble(all_mut_A)
#T nucleotide (changing to U for RNA)
all_mut_U <-subset(all_mut, RefBase =="T")
all_mut_U <- as_tibble(all_mut_U)
#C nucleotide
all_mut_C <-subset(all_mut, RefBase =="C")
all_mut_C <- as_tibble(all_mut_C)

## add columns for each base change 
all_mut_G <- all_mut_G %>% mutate(
  G_to_A = (Acnt / Coverage),
  G_to_C = (Ccnt / Coverage),
  G_to_U = (Tcnt / Coverage)
)

all_mut_A <- all_mut_A %>% mutate(
  A_to_G = (Gcnt / Coverage),
  A_to_C = (Ccnt / Coverage),
  A_to_U = (Tcnt / Coverage)
)

all_mut_U <- all_mut_U %>% mutate(
  U_to_G = (Gcnt / Coverage),
  U_to_A = (Acnt / Coverage),
  U_to_C = (Ccnt / Coverage)
)

all_mut_C <- all_mut_C %>%mutate(
  C_to_G = (Gcnt / Coverage),
  C_to_A = (Acnt / Coverage),
  C_to_U = (Tcnt / Coverage)
)

#melts mutation rows to be variables in "basechange" with their value next to it
all_mut_A <-all_mut_A %>% gather(basechange, value, "A_to_G":"A_to_U")
all_mut_C <-all_mut_C %>% gather(basechange, value, "C_to_G":"C_to_U")
all_mut_U <-all_mut_U %>% gather(basechange, value, "U_to_G":"U_to_C")
all_mut_G <-all_mut_G %>% gather(basechange, value, "G_to_A":"G_to_U")

basechange_t <-bind_rows(all_mut_A, all_mut_C, all_mut_G, all_mut_U)
#filter only positions that have min 200 coverage
basechange_t <- subset(basechange_t, Coverage > 200)


#getting mean  basechange value per sample 
meanGU <- basechange_t %>% filter(basechange=="G_to_U") %>% group_by(Cohort, DPI, Sample.Type, Animal.label) %>% 
  summarize(G_to_U=mean(value, na.rm=T ), stdev=sd(value, na.rm=T))
meanGC <- basechange_t %>% filter(basechange=="G_to_C") %>% group_by(Cohort, DPI, Sample.Type, Animal.label) %>% 
  summarize(G_to_C=mean(value, na.rm=T ), stdev=sd(value, na.rm=T))
meanGA <- basechange_t %>% filter(basechange=="G_to_A") %>% group_by(Cohort, DPI, Sample.Type, Animal.label) %>% 
  summarize(G_to_A=mean(value, na.rm=T ), stdev=sd(value, na.rm=T)) 

meanCU <- basechange_t %>% filter(basechange=="C_to_U") %>% group_by(Cohort, DPI, Sample.Type, Animal.label) %>% 
  summarize(C_to_U=mean(value, na.rm=T ), stdev=sd(value, na.rm=T)) 
meanCA <- basechange_t %>% filter(basechange=="C_to_A") %>% group_by(Cohort, DPI, Sample.Type, Animal.label) %>% 
  summarize(C_to_A=mean(value, na.rm=T ), stdev=sd(value, na.rm=T)) 
meanCG <- basechange_t %>% filter(basechange=="C_to_G") %>% group_by(Cohort, DPI, Sample.Type, Animal.label) %>% 
  summarize(C_to_G=mean(value, na.rm=T ), stdev=sd(value, na.rm=T)) 

meanAU <- basechange_t %>% filter(basechange=="A_to_U") %>% group_by(Cohort, DPI, Sample.Type, Animal.label) %>% 
  summarize(A_to_U=mean(value, na.rm=T ), stdev=sd(value, na.rm=T)) 
meanAC <- basechange_t %>% filter(basechange=="A_to_C") %>% group_by(Cohort, DPI, Sample.Type, Animal.label) %>% 
  summarize(A_to_C=mean(value, na.rm=T ), stdev=sd(value, na.rm=T)) 
meanAG <- basechange_t %>% filter(basechange=="A_to_G") %>% group_by(Cohort, DPI, Sample.Type, Animal.label) %>% 
  summarize(A_to_G=mean(value, na.rm=T ), stdev=sd(value, na.rm=T)) 

meanUC <- basechange_t %>% filter(basechange=="U_to_C") %>% group_by(Cohort, DPI, Sample.Type, Animal.label) %>% 
  summarize(U_to_C=mean(value, na.rm=T ), stdev=sd(value, na.rm=T)) 
meanUG <- basechange_t %>% filter(basechange=="U_to_G") %>% group_by(Cohort, DPI, Sample.Type, Animal.label) %>% 
  summarize(U_to_G=mean(value, na.rm=T ), stdev=sd(value, na.rm=T)) 
meanUA <- basechange_t %>% filter(basechange=="U_to_A") %>% group_by(Cohort, DPI, Sample.Type, Animal.label) %>% 
  summarize(U_to_A=mean(value, na.rm=T ), stdev=sd(value, na.rm=T)) 

#put back together for plotting 
mean_GCGU <- merge(meanGC, meanGU,
                   by=c("Cohort", "DPI", "Sample.Type", "Animal.label"))
mean_CACG <- merge(meanCA, meanCG,
                   by=c("Cohort", "DPI", "Sample.Type", "Animal.label"))
mean_UCAU <- merge(meanUC, meanAU,
                   by=c("Cohort", "DPI", "Sample.Type", "Animal.label"))
mean_AGAC <- merge(meanAG, meanAC,
                   by=c("Cohort", "DPI", "Sample.Type", "Animal.label"))
mean_UAUG <- merge(meanUA, meanUG,
                   by=c("Cohort", "DPI", "Sample.Type", "Animal.label"))

meanGACU <- merge(meanGA, meanCU,
                  by=c("Cohort", "DPI", "Sample.Type", "Animal.label"))

mean_basechange <- merge(merge(merge(merge(merge(
  mean_GCGU, 
  mean_CACG, by=c("Cohort", "DPI", "Sample.Type", "Animal.label")),
  mean_UCAU, by=c("Cohort", "DPI", "Sample.Type", "Animal.label")),
  mean_AGAC, by=c("Cohort", "DPI", "Sample.Type", "Animal.label")),
  mean_UAUG, by=c("Cohort", "DPI", "Sample.Type", "Animal.label")),
  meanGACU, by=c("Cohort", "DPI", "Sample.Type", "Animal.label"))

melt_allbc <- reshape2::melt(mean_basechange, na.rm=F, id.vars=c("Cohort", "DPI", "Sample.Type", "Animal.label"), 
                             measure.vars=c("G_to_A", "G_to_U", "G_to_C", "C_to_G", "C_to_A", "C_to_U", 
                                            "A_to_G", "A_to_U", "A_to_C", "U_to_A", "U_to_G", "U_to_C"))

melt_allbc$variable <- factor(melt_allbc$variable)
melt_allbc$basechange <- factor(melt_allbc$basechange, levels=c("G_to_A", "G_to_U", "G_to_C", "C_to_U", "C_to_A", "C_to_G", 
                                                                "A_to_G", "A_to_U", "A_to_C", "U_to_A", "U_to_G", "U_to_C"))

colnames(melt_allbc)[which(names(melt_allbc) == "variable")] <- "basechange"

stats4 <- melt_allbc %>%
  filter(Sample.Type=="Swab") %>% 
  filter(DPI != "5") %>% 
  filter(DPI != "4") %>%
  rstatix::group_by(Cohort) %>% 
  rstatix::wilcox_test(value~DPI) %>% 
  adjust_pvalue(method="none") %>% 
  add_significance("p") %>%
  add_xy_position(x = "Cohort", dodge=0.8) %>%
  remove_ns(col="p")

base.change.all.plot <- ggplot(melt_allbc %>% filter(Sample.Type == "Swab") %>% 
                                 filter(DPI != "5") %>% filter(DPI != "4"), 
                               aes(x=Cohort, y=value, colour= as.factor(DPI))) +
  geom_boxplot(outlier.size = 0.1) +
  ylab("Mean proportion of base change") + xlab("") +
  ggtitle("All base changes") +
  theme(axis.text.x = element_text( vjust = 0.5, hjust=1))+
  facet_wrap(~basechange, ncol=3) + 
  scale_y_log10() + ylim(0,0.003) +
  theme_bw() +
  scale_colour_viridis_d(direction=-1, name="DPI")+
  theme(axis.text.y = element_text(colour="black", size=12),
        axis.title.y = element_text(size = 12, face="bold"),
        axis.text.x = element_text(size=12, colour="black", angle=45, hjust=1),
        axis.title.x = element_text(size = 12, face="bold"),
        strip.background = element_rect(size=0.5, linetype="solid", fill="white"),
        strip.text.x = element_text(size = 12, color = "black", face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))#+

base.change.all.plot

#ggsave(base.change.all.plot, filename = "~/projects/agile_mpv_seq_project/Reports/Agile_paper_figs/sup1_allbase.png", device="png", dpi=300, width=8, height=12)

#isolating GACU for main plots
meanGACU <- merge(meanGA, meanCU, by =c("Cohort", "DPI", "Sample.Type", "Animal.label"))
melt_meanGACU <- reshape2::melt(meanGACU, na.rm=F, id.vars=c("Cohort", "DPI", "Sample.Type", "Animal.label"), measure.vars=c("G_to_A", "C_to_U"))
melt_meanGACU$variable <- factor(melt_meanGACU$variable, levels=c("G_to_A", "C_to_U"))
melt_meanGACU <- rename(melt_meanGACU, basechange=variable)


#plotting 
#setting labels for facet 
base.labs <- c("C to U", " G to A")
names(base.labs) <- c("C_to_U", "G_to_A")

  
colnames(melt_meanGACU)[which(names(melt_meanGACU) == "variable")] <- "basechange"

stats3 <- melt_meanGACU %>%
    filter(Sample.Type=="Swab") %>% 
    rstatix::group_by(Cohort,basechange) %>% 
    rstatix::wilcox_test(value~DPI) %>% 
    adjust_pvalue(method="none") %>% 
    add_significance("p") %>%
    add_xy_position(x = "DPI", dodge=0.8) %>%
    remove_ns(col="p")

base.change.plot <- ggplot(melt_meanGACU %>% filter(Sample.Type == "Swab") %>% 
                             filter(DPI != "5") %>% filter(DPI != "4"), 
                           aes(x=as.factor(DPI), y=value, colour= as.factor(DPI))) +
  geom_boxplot(outlier.size = 0.1, varwidth = TRUE) +
  ylab("Mean proportion of base change") + xlab("Day of swab sample") +
  facet_grid(basechange~Cohort, labeller = labeller(basechange = base.labs, Cohort = label_wrap_gen(10))) +
  theme_bw() +
  ylim(0,0.0026)+
  scale_colour_viridis_d(direction=1, name="DPI")+
  theme(axis.text.y = element_text(colour="black", size=12),
        axis.title.y = element_text(size = 12, face="bold"),
        axis.text.x = element_text(size=12, colour="black", angle=45, hjust=1),
        axis.title.x = element_text(size = 12, face="bold"),
        strip.background = element_rect(size=0.5, linetype="solid", fill="white"),
        strip.text.x = element_text(size = 12, color = "black", face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  stat_pvalue_manual(stats3, label = "p.signif", y.position = 0.0025, tip.length = 0)
base.change.plot 



#png for git 
ggsave(base.change.plot, filename="~/projects/cyclophosphamide-mice/nimagen2/figures2/basechange-cu-ga-plot-stats.png", device="png", dpi=300, width=13, height=6)
#tiff for pub
ggsave(base.change.plot, filename="~/projects/cyclophosphamide-mice/nimagen2/figures2/basechange-cu-ga-plot-stats.tiff", device="tiff", dpi=300, width=13, height=6)
