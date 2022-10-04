# Script used to plot the S: H665Y mutation in throat swabs and tissues derived from mice with and without chemical immunosupression
# with and without medical countermeasure Molnupiravir and PF-07321332 (Pfizer compound)

library(tidyverse)
library(ggplot2)

setwd("~/projects/cyclophosphamide-mice/nimagen2/DiversiTools/")


##### Read in diversitools *_AA.txt output
diversitools_readin_aa <- function(filepath ="./", file_pattern, file_list) {
  
  temp <- list.files(filepath, pattern=file_pattern) # creates list with filenames as in the directory
  
  file_list <- list() # creates empty list to add to in the for loop
  
  for (i in 1:length(temp)) { # read in all data files as DFs, add Proportion(nonsyn) column to all 
    
    file_list[[i]] <- read.delim(temp[i]) 
    file_list[[i]]$AAcoverage <- as.numeric(file_list[[i]]$AAcoverage)
    file_list[[i]]$TopAAcnt <- as.numeric(file_list[[i]]$TopAAcnt) # changes class to numeric
    file_list[[i]]$SndAAcnt <- as.numeric(file_list[[i]]$SndAAcnt)
    file_list[[i]]$TrdAAcnt <- as.numeric(file_list[[i]]$TrdAAcnt)
    file_list[[i]]$ProportionTopAA <-  file_list[[i]]$TopAAcnt / file_list[[i]]$AAcoverage
    file_list[[i]]$ProportionSndAA <-  file_list[[i]]$SndAAcnt / file_list[[i]]$AAcoverage
    file_list[[i]]$ProportionTrdAA <-  file_list[[i]]$TrdAAcnt / file_list[[i]]$AAcoverage
    file_list[[i]]$AAPosition <- as.numeric(file_list[[i]]$AAPosition)
    file_list[[i]][,7:14] <- lapply(file_list[[i]][,7:14], as.character)
    file_list[[i]]$TopCodoncnt <- as.numeric(file_list[[i]]$TopCodoncnt) # changes class to numeric
    file_list[[i]]$SndCodoncnt <- as.numeric(file_list[[i]]$SndCodoncnt)
    file_list[[i]]$TrdCodoncnt <- as.numeric(file_list[[i]]$TrdCodoncnt)
    file_list[[i]]$AAPsn <- 1:nrow(file_list[[i]]) #adding another column to set the positions as the number of rows so that they are not filtered out later on
    file_list[[i]]$Protein <- factor(file_list[[i]]$Protein, levels=c("nsp1", "nsp2", "nsp3", "nsp4", "nsp5", "nsp6", "nsp7", "nsp8", "nsp9",
                                                                      "nsp10", "nsp12", "nsp12_2", "nsp13", "nsp14", "nsp15", "nsp16", "S", 
                                                                      "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10"))
    #adding levels to the Protein column so that they come up in the right order along the output graph and in the key
    file_list[[i]] <- file_list[[i]]  %>% filter(Sample != "Sample") %>% filter(AAcoverage > 20) #%>% filter(RefAA != TopAA)
    cols <- c("Protein", "RefAA", "AAPosition", "TopAA") # creates vector of cols we want to paste in next step
    file_list[[i]]$subsTopAA <- do.call(paste, c(file_list[[i]][cols])) # creates new column "mutant" pasting the protein name, refAA, AAPosition, TopAA
    # #write.csv(filtered_AA, file = paste0(outdir, names(myfilenames)[i], "_amino.acids.test.csv")) # write out as csv into your chosen directory
    # 
    cols2 <- c("Protein", "RefAA", "AAPosition", "SndAA")
    file_list[[i]]$subsSndAA <- do.call(paste, c(file_list[[i]][cols2]))
    cols3 <- c("Protein", "RefAA", "AAPosition", "TrdAA")
    file_list[[i]]$subsTrdAA <- do.call(paste, c(file_list[[i]][cols3]))
    }
  
  myfilenames <- gsub("_AA.txt","", temp) # give names to each file, removing unnecessary strings
  names(file_list) <- myfilenames
  return(file_list)
}


aa_list_sp <- diversitools_readin_aa(filepath = "./", file_pattern = "_AA.txt", file_list = aa_list)

for (i in 1:length(aa_list_sp)) {
  aa_list_sp[[i]]$Sample <- gsub(".final.bam", "", aa_list_sp[[i]]$Sample)
}

all_aa_sp <- bind_rows(aa_list_sp)

#read in metadata
metadata_sp <- read.csv("../nextclade/metadata3.csv")

#change colname to match dataframe
colnames(metadata_sp)[which(names(metadata_sp) == "seqName")] <- "Sample"

all_aa_sp <- merge(all_aa_sp, metadata_sp, by="Sample")

all_aa_sp <- all_aa_sp %>% filter(Exclude!="Y")

#code to check class of entries and identify unique values - highlights any spelling mistakes or unexpected groups
sapply(melt_all_aa_sp, class)
unique(melt_all_aa_sp$DPI)
unique(melt_all_aa_sp$Sample.Type)

#reshape the data

melt_all_aa_sp <- reshape2::melt(all_aa_sp, id.vars = c('AAPosition', "AAPsn", "Sample", "clade", "Animal.label", "DPI", "Sample.Type", "Cohort", "Protein", "TopAA", "SndAA", "TrdAA" ),
                                 measure.vars=  c("ProportionTopAA", "ProportionSndAA", "ProportionTrdAA"))
melt_all_aa_sp$variable <- as.factor(melt_all_aa_sp$variable)
melt_all_aa_sp$DPI <- factor(melt_all_aa_sp$DPI, levels=c("1", "3", "4", "6", "7"))
melt_all_aa_sp$Sample.Type <- factor(melt_all_aa_sp$Sample.Type, levels = c("Swab", "Lung RNA","Nasal Tissue"))
melt_all_aa_sp$Cohort <- factor(melt_all_aa_sp$Cohort, levels=c("Vehicle", "Cyclophosphamide only", "Molnupiravir","Cyclophosphamide and Molnupiravir","Cyclophosphamide and Pfizer","Cyclophosphamide, Pfizer and Molnupiravir"))

levels(melt_all_aa_sp$Cohort)

#just S 655 labels ####
all_aa_sp_655 <- all_aa_sp %>% filter(Protein=="S") %>% filter(AAPosition =="655")
mut1_AA<- all_aa_sp_655[,"TopAA"]
mut2_AA<- all_aa_sp_655[, "SndAA"]
mut3_AA<- all_aa_sp_655[, "TrdAA"]
mut_label_AA_655 <- c(mut1_AA, mut2_AA, mut3_AA)

#### summarise for S:655 
S_655 <- select(all_aa_sp, 3:4, 14, 16, 18, 27:29, 34:42) %>% filter(Protein=="S") %>% filter(AAPosition =="655")
#melt AA first
melt_S655_ <- pivot_longer(S_655, cols=c("TopAA", "SndAA", "TrdAA"), names_to = "AA_rank", values_to = "AA")
#then melt proportions
melt_S655_1 <- pivot_longer(melt_S655_, cols=c("ProportionTopAA", "ProportionSndAA", "ProportionTrdAA"), names_to = "Prop_AA_name", values_to = "Prop_AA")
#melting props gives a 3 rows per AA so changing Prop_AA_names so that they match format of AA_rank
melt_S655_1$Prop_AA_name <- gsub("Proportion", "", melt_S655_1$Prop_AA_name)
#keep only rows where AA_rank and Prop_AA_name match to get only the value that matches the rank 
melt_S655_final <- melt_S655_1 %>% filter(AA_rank == Prop_AA_name)
melt_S655_final$Cohort <- factor(melt_S655_final$Cohort, levels=c("Vehicle", "Cyclophosphamide only", "Molnupiravir","Cyclophosphamide and Molnupiravir","Cyclophosphamide and Pfizer","Cyclophosphamide, Pfizer and Molnupiravir"))


plot_S655_mv_line <- ggplot(melt_S655_final %>% filter(AA_rank!="TrdAA") %>% filter(AA!="<NA>"), aes(x=DPI, y=Prop_AA, colour=AA, shape = Animal.label)) + 
  geom_point(size=2.5) + geom_line() +
  facet_grid(Sample.Type~Cohort, 
             labeller = labeller(Cohort = label_wrap_gen(10))) +
  #scale_x_continuous(breaks=c(0,1,4,6,7), limits = c(0,7)) +
  ylim(0,1) + ylab("Proportion of amino acid at S:655\n") + xlab("\nDPI") +
  scale_colour_viridis_d(begin = 0, end = 0.8, direction=-1) +
  theme_bw() +
  labs(colour="Amino\nacid") +
  labs(shape="Animal") +
  theme(axis.text = element_text(colour="black", size=12, vjust=0.2),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold", vjust=2),
        strip.background = element_rect(size=0.5, linetype="solid", fill="white"),
        strip.text = element_text(size = 12, color = "black", face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(face="bold", hjust=0.5),
        legend.text = element_text(size=10))

plot_S655_mv_line

#write out
#png for git 
ggsave(plot_S655_mv_line, filename=("~/projects/cyclophosphamide-mice/nimagen2/figures2/S655_mv_line.png"), device="png", dpi=300, width=14, height=6)
#tiff for pub 
ggsave(plot_S655_mv_line, filename=("~/projects/cyclophosphamide-mice/nimagen2/figures2/S655_mv_line.tiff"), device="tiff", dpi=300, width=14, height=6)


#### just S 655 labels ####

all_aa_sp_orf8 <- all_aa_sp %>% filter(Protein=="ORF8") %>% filter(AAPosition =="27")
mut1_AA<- all_aa_sp_orf8[,"TopAA"]
mut2_AA<- all_aa_sp_orf8[, "SndAA"]
mut3_AA<- all_aa_sp_orf8[, "TrdAA"]
mut_label_AA_655 <- c(mut1_AA, mut2_AA, mut3_AA)

#summarise for S:655
S_655 <- select(all_aa_sp, 3:4, 14, 16, 18, 27:29, 34:42) %>% filter(Protein=="ORF8") %>% filter(AAPosition =="27")
#melt AA first
melt_S655_ <- pivot_longer(S_655, cols=c("TopAA", "SndAA", "TrdAA"), names_to = "AA_rank", values_to = "AA")
#then melt proportions
melt_S655_1 <- pivot_longer(melt_S655_, cols=c("ProportionTopAA", "ProportionSndAA", "ProportionTrdAA"), names_to = "Prop_AA_name", values_to = "Prop_AA")
#melting props gives a 3 rows per AA so changing Prop_AA_names so that they match format of AA_rank
melt_S655_1$Prop_AA_name <- gsub("Proportion", "", melt_S655_1$Prop_AA_name)
#keep only rows where AA_rank and Prop_AA_name match to get only the value that matches the rank 
melt_S655_final <- melt_S655_1 %>% filter(AA_rank == Prop_AA_name)



plot_S655_mv_line <- ggplot(melt_S655_final %>% filter(AA_rank!="TrdAA") %>% filter(AA!="<NA>"), aes(x=DPI, y=Prop_AA, colour=AA, shape = Animal.label)) + 
  geom_point(size=2.5) + geom_line() +
  facet_grid(Sample.Type~Cohort, 
             labeller = labeller(Cohort = label_wrap_gen(10))) +
  ylim(0,1) + ylab("Proportion of amino acid at S:655\n") + xlab("\nDPI") +
  scale_colour_viridis_d(begin = 0, end = 0.8, direction=-1) +
  theme_bw() +
  labs(colour="Amino\nacid") +
  labs(shape="Animal") +
  theme(axis.text = element_text(colour="black", size=12, vjust=0.2),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold", vjust=2),
        strip.background = element_rect(size=0.5, linetype="solid", fill="white"),
        strip.text = element_text(size = 12, color = "black", face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(face="bold", hjust=0.5),
        legend.text = element_text(size=10))

plot_S655_mv_line

#write out
#png for git 
ggsave(plot_S655_mv_line, filename=("~/projects/cyclophosphamide-mice/nimagen/figures/S655_mv_line.png"), device="png", dpi=300, width=14, height=6)
#tiff for pub 
ggsave(plot_S655_mv_line, filename=("~/projects/cyclophosphamide-mice/nimagen/figures/S655_mv_line.tiff"), device="tiff", dpi=300, width=14, height=6)




