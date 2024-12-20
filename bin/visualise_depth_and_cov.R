# visualise coverage and depth
library(dplyr)

data<-read.csv("../bam/coverage_results.csv")
metadata<-read.csv("../../../viral-genomics-immunosupression-and-countermeasures/metadata/metadata3.csv")%>%
  filter(Exclude!="Y")

metadata$Cohort <- gsub("Pfizer", "Nirmatrelvir",metadata$Cohort)

metadata$Filename <- paste(metadata$seqName,".final.bam")
metadata$Filename <- gsub(" ", "", metadata$Filename)

metadata$Filename

colnames(metadata)

data<-data %>% filter(Filename %in% metadata$Filename)

metadata_clean <- metadata %>%
  select(c(Filename, Cohort, Animal, Sample.Type,
           DPI,ct,Animal,seqName))

data_merged<-left_join(metadata_clean, data, by="Filename")

data_merged$unique_sample <- paste(data_merged$Cohort, data_merged$Animal, 
                                   "DPI:",data_merged$DPI)

colnames(data_merged)

library(ggplot2)

swab<-data_merged%>%filter(Sample.Type=="Swab")%>%
  filter(Threshold==10)

min<-min(swab$Average_Depth)
max<-max(swab$Average_Depth)

# Define custom Cohort order
custom_cohort_order <- c(
  "Vehicle", 
  "Cyclophosphamide only", 
  "Molnupiravir", 
  "Cyclophosphamide and Molnupiravir",
  "Cyclophosphamide and Nirmatrelvir",
  "Cyclophosphamide, Nirmatrelvir and Molnupiravir"
)

# Preprocess the data
data_processed <- data_merged %>%
  filter(Sample.Type == "Swab", Threshold == 100) %>%
  mutate(Cohort = factor(Cohort, 
                         levels = custom_cohort_order)) %>% # Set custom Cohort order
  group_by(Cohort) %>% # Group by Cohort
  arrange(Cohort, DPI, Animal) %>% # Sort within each Cohort by DPI and animal_id
  ungroup() %>%
  mutate(unique_sample = factor(unique_sample, levels = unique(unique_sample)))%>%
  mutate(unique_sample = fct_rev(unique_sample))# Reorder y-axis

#plot average depth of each sample
depth<-ggplot(data_processed %>%
         filter(Sample.Type=="Swab")%>%
         filter(Threshold==100), aes(Average_Depth,unique_sample))+
  geom_bar(stat="identity", fill="#1A85FF")+
  theme_classic()+
  xlab("log10(average depth)")+
  geom_vline(xintercept = min, linetype="dashed",colour="#D41159")+
  geom_vline(xintercept= max,linetype="dashed",colour="#D41159")+
  scale_x_log10(labels=scales::comma)+
  annotate("label",label= "min: 1845.89", x= 1000, 
           y = 25, angle = 90, size=3,colour="#D41159")+
  annotate("label",label= "max: 57027.32", x= 100000, 
           y = 25, angle = 90, size=3, colour="#D41159")+
  ylab("")

ggsave("../figures2/depth.tiff", 
       device = "tiff",
       width = 8,
       height = 7,
       units = "in",
       dpi=300)

cov<-ggplot(data_processed%>%
         filter(Sample.Type=="Swab")%>%
         filter(Threshold==100),aes(Breadth_of_Coverage,unique_sample))+
  geom_bar(stat="identity",fill="#1A85FF")+
  geom_vline(xintercept = 90, linetype = "dashed", colour="#D41159")+
  xlab("% coverage 100X")+
  ylab("")+
  theme_classic()

ggsave("../figures2/depth_100x.tiff", device = "tiff",
       width = 8,
       height = 7,
       units = "in",
       dpi=300)

plot_grid(cov + theme(plot.margin = unit(c(1, 1, 1, 1), "lines")),
          depth + theme(axis.ticks.y=element_blank(),
                      axis.text.y = element_blank(),
                      plot.margin = unit(c(1, 1, 1, 1), "lines")),
          labels = "AUTO",
          rel_widths = c(1,0.5),
          align = "h")


ggsave("../figures2/depth_cov.tiff", 
       device = "tiff",
       width = 10,
       height = 7,
       units = "in",
       dpi=300)
