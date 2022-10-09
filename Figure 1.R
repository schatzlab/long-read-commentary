# plot sequencing stats
# Shujun Ou (shujun.ou.1@gmail.com)
# 10/01/2020

library(tidyr)
library(ggplot2)
library(wesanderson)
library(scales)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(patchwork)
require(data.table)
library(stringr)

workdir = '/Users/oushujun/My Drive/study/JHU/NatMethods commentary/'
setwd(workdir)

## color settings
vendor_colors = c("PacBio"="#d7228c", "Nanopore"="#0a5d74", "Illumina"="#f6932a", 'Others'='black')

# read and process data
euk_genomes = read.table('ncbi_eukaryota_genomes.tsv', sep = '\t', fill = T, header = T)
euk_genomes$Assembly.Submission.Date = as.Date(euk_genomes$Assembly.Submission.Date)
euk_genomes$Assembly.Submission.Year = sub('-.*', '', euk_genomes$Assembly.Submission.Date)
seq_tech = read.csv('sequencing_tech.csv', header = T)
seq_tech$Year = sub('.*/', '20', seq_tech$Year)
seq_tech$Mean.read.length = sub(',', '', seq_tech$Mean.read.length)
seq_tech$Mean.read.length = as.numeric(seq_tech$Mean.read.length)
seq_tech$Vendor[!str_detect('Illumina|PacBio|Nanopore', seq_tech$Vendor)] = 'Ohters'
seq_tech$Vendor = as.factor(seq_tech$Vendor)
seq_tech$Year = as.numeric(seq_tech$Year)
ont_max = read.csv('Nanopore_max_len.csv')
ont_max$Length = as.numeric(sub(',', '', ont_max$Length))
ont_max$Date = as.numeric(sub('.*\\/', '', ont_max$Date))
ont_max = ont_max %>% group_by(Date) %>% summarise(Year = Date, max.length = max(Length, na.rm=T))

# plot sequencing technologies
seq_tech_length_plot = ggplot(subset(seq_tech, Year >= 2002), 
                              aes(Year, Mean.read.length/1000, color = Vendor)) + 
  scale_color_manual(values = vendor_colors) + xlim(c(2001, 2022)) +
  geom_point() + labs(y = 'Mean Read Length (kbp)') +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = 'grey', size = 0.2, linetype = 1),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = c(0.5, 0.7))
seq_tech_length_plot

# plot max ont length
ont_max_plot = ggplot(ont_max, aes(Year, max.length/1000000)) + geom_line(color = '#0a5d74') +
  labs(y = 'Longest Read Length (Mb)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background =element_rect(colour ="black", size=0.5))

# make insert
seq_combined_plot = seq_tech_length_plot + 
  inset_element(ont_max_plot, left = 0.05, bottom = 0.25, right = 0.4, top = 0.98)
seq_combined_plot


# find most popular species
head(sort(table(euk_genomes$Organism.Name), decreasing = T), 100)
head(sort(table(subset(euk_genomes, Assembly.Stats.Total.Sequence.Length > 2000000000)$Organism.Name), decreasing = T), 100)
target_list = c('Homo sapiens', 'Mus musculus', 'Ovis aries', 'Zea mays', 'Triticum aestivum', 'Hordeum vulgare')

# merge species names
subset(euk_genomes, str_detect(euk_genomes$Organism.Name, 'Hordeum vulgare'))$Organism.Name
subset(euk_genomes, str_detect(euk_genomes$Organism.Name, 'Hordeum vulgare'))
euk_genomes$Organism.Name[euk_genomes$Organism.Name=='Zea mays subsp. mays'] = 'Zea mays'
euk_genomes$Organism.Name[euk_genomes$Organism.Name=='Mus musculus musculus'] = 'Mus musculus'
euk_genomes$Organism.Name[euk_genomes$Organism.Name=='Hordeum vulgare subsp. vulgare'] = 'Hordeum vulgare'

# filter out abnormal data
subset(euk_genomes, Organism.Name == 'Homo sapiens' & Assembly.Submission.Year < 2005)
euk_genomes = subset(euk_genomes, !is.na(Assembly.Submission.Year))
euk_genomes = subset(euk_genomes, Assembly.Stats.Contig.N50 < 500000000)
euk_genomes = subset(euk_genomes, Assembly.Stats.Total.Sequence.Length > 1000000)
euk_genomes = subset(euk_genomes, !(Assembly.Stats.Total.Sequence.Length < 2000000000 &
                                      Organism.Name == 'Homo sapiens'))
euk_genomes = subset(euk_genomes, !(Assembly.Stats.Total.Sequence.Length > 5000000000 &
                                      Organism.Name == 'Homo sapiens'))
euk_genomes = subset(euk_genomes, !(Assembly.Stats.Total.Sequence.Length < 2000000000 &
                                      Organism.Name == 'Zea mays'))
euk_genomes = subset(euk_genomes, !(Assembly.Stats.Total.Sequence.Length < 2000000000 &
                                      Organism.Name == 'Triticum aestivum'))
euk_genomes = subset(euk_genomes, !(Assembly.Stats.Total.Sequence.Length < 2000000000 &
                                      Organism.Name == 'Hordeum vulgare'))

#mean N50 ~ year
options("scipen"=100, "digits"=1)
euk_genomes_mean = subset(euk_genomes, Organism.Name %in% target_list) %>% 
  group_by(Organism.Name, Assembly.Submission.Year, .add = TRUE) %>%
  summarise(contig.N50_mean = mean(Assembly.Stats.Contig.N50, na.rm=T))
euk_genomes_mean$Assembly.Submission.Year = as.numeric(euk_genomes_mean$Assembly.Submission.Year)

# add genome size to genome names
new_name = data.frame(species = factor(c('Triticum aestivum', 'Hordeum vulgare', 'Zea mays',
                                     'Homo sapiens', 'Mus musculus', 'Ovis aries')),
                   species_new = factor(c('Triticum aestivum (~16 Gb)', 'Hordeum vulgare (~5.3 Gb)', 'Zea mays (2.4 Gb)',
                                          'Homo sapiens (3.05 Gb)', 'Mus musculus (2.6 Gb)', 'Ovis aries (~2.6 Gb)')))
euk_genomes_mean$Organism.Name.size = new_name$species_new[match(euk_genomes_mean$Organism.Name, new_name$species)]
euk_genomes_mean$Organism.Name.size = factor(euk_genomes_mean$Organism.Name.size, 
                                        levels = c('Triticum aestivum (~16 Gb)', 'Hordeum vulgare (~5.3 Gb)', 'Zea mays (2.4 Gb)',
                                                   'Homo sapiens (3.05 Gb)', 'Mus musculus (2.6 Gb)', 'Ovis aries (~2.6 Gb)'))

# function to increase vertical spacing between legend keys
# @clauswilke
draw_key_polygon3 <- function(data, params, size) {
  lwd <- min(data$size, min(size) / 4)
  
  grid::rectGrob(
    width = grid::unit(0.6, "npc"),
    height = grid::unit(0.6, "npc"),
    gp = grid::gpar(
      col = data$colour,
      fill = alpha(data$fill, data$alpha),
      lty = data$linetype,
      lwd = lwd * .pt,
      linejoin = "mitre"
    ))
}

# make a plot
euk_genomes_mean_plot = ggplot(euk_genomes_mean, aes(Assembly.Submission.Year, color = Organism.Name.size, 
                                                     contig.N50_mean/1000000, shape = Organism.Name.size,)) + 
  geom_point() +
  geom_line() +
  scale_y_log10() + xlim(c(2001, 2022)) +
  labs(x = 'Assembly Submission Year', y = 'Mean Contig N50 (Mb)', color = 'Species', shape = 'Species') +
  geom_vline(xintercept = 2001.7, color = 'black', linetype = 'dotted') +
  geom_vline(xintercept = 2007, color = '#f6932a', linetype = 'dotted') +
  geom_vline(xintercept = 2010, color = '#d7228c', linetype = 'dotted') +
  geom_vline(xintercept = 2014, color = '#0a5d74', linetype = 'dotted') +
  geom_vline(xintercept = 2015, color = '#d7228c', linetype = 'dotted') +
  geom_vline(xintercept = 2018, color = '#0a5d74', linetype = 'dotted') +
  geom_vline(xintercept = 2019, color = '#d7228c', linetype = 'dotted') +
  geom_vline(xintercept = 2021, color = '#0a5d74', linetype = 'dotted') +
  annotate('text', x = 2001.4, y = 100, label = 'Sanger/454', angle = 90, hjust = 1, size = 3) +
  annotate('text', x = 2006.7, y = 100, label = 'Illumina', angle = 90, hjust = 1, size = 3) +
  annotate('text', x = 2009.7, y = 100, label = 'PB RS', angle = 90, hjust = 1, size = 3) +
  annotate('text', x = 2013.7, y = 100, label = 'ONT MinION', angle = 90, hjust = 1, size = 3) +
  annotate('text', x = 2014.7, y = 100, label = 'PB Sequel', angle = 90, hjust = 1, size = 3) +
  annotate('text', x = 2017.7, y = 1e-3, label = 'ONT PromethION', angle = 90, hjust = 0, size = 3) +
  annotate('text', x = 2018.7, y = 1e-3, label = 'PB HiFi', angle = 90, hjust = 0, size = 3) +
  annotate('text', x = 2020.7, y = 1e-3, label = 'ONT Duplex', angle = 90, hjust = 0, size = 3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = c(0.18, 0.78), legend.key = element_rect(fill = NA), legend.title = element_blank(), 
        legend.background = element_rect(fill = NA), legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, 'cm'))
euk_genomes_mean_plot


# combine plots
pdf("Figure 1.pdf", width=10.5,height=7.5,pointsize=12, paper='special')
seq_combined_plot / euk_genomes_mean_plot + 
  plot_annotation(tag_levels = 'a')  & theme(plot.tag = element_text(size = 16, face = 'bold'))
dev.off()



##########################
### Other unused plots ###
##########################
euk_genome_contigN50_plot_ori = ggplot(euk_genomes) + 
  geom_boxplot(aes(Assembly.Submission.Year, Assembly.Stats.Contig.N50/1000000)) +
  annotate('text', x = 6, y = 250, label = 'SOLiD') +
  annotate('text', x = 13, y = 250, label = 'P5C3') +
  annotate('text', x = 18, y = 250, label = 'R9.4') +
  annotate('text', x = 19, y = 250, label = 'CCS') +
  labs(title = 'Eukaryota genomes', y = 'Contig N50 (Mb)') +
  theme_bw()

euk_genome_contigN50_plot_log = ggplot(euk_genomes) + 
  geom_boxplot(aes(Assembly.Submission.Year, Assembly.Stats.Contig.N50/1000000)) +
  annotate('text', x = 6, y = 250, label = 'SOLiD') +
  annotate('text', x = 13, y = 250, label = 'P5C3') +
  annotate('text', x = 18, y = 250, label = 'R9.4') +
  annotate('text', x = 19, y = 250, label = 'CCS') +
  scale_y_log10() + 
  labs(title = 'Eukaryota genomes', y = 'Contig N50 (Mb)') +
  theme_bw()

euk_genome_contigN50_plot_ori / euk_genome_contigN50_plot_log


human_genome_contigN50_plot_ori = ggplot(subset(euk_genomes, Organism.Name == 'Homo sapiens')) + 
  geom_boxplot(aes(Assembly.Submission.Year, Assembly.Stats.Contig.N50/1000000)) +
  labs(title = 'Homo sapiens', y = 'Contig N50 (Mb)') +
  theme_bw()

human_genome_contigN50_plot_log = ggplot(subset(euk_genomes, Organism.Name == 'Homo sapiens')) + 
  geom_boxplot(aes(Assembly.Submission.Year, Assembly.Stats.Contig.N50/1000000)) +
  scale_y_log10() + 
  labs(title = 'Homo sapiens', y = 'Contig N50 (Mb)') +
  theme_bw()

human_genome_contigN50_plot_ori / human_genome_contigN50_plot_log


euk_genome_scfN50_plot = ggplot(euk_genomes) + 
  geom_boxplot(aes(Assembly.Submission.Year, Assembly.Stats.Scaffold.N50/1000000)) +
  labs(title = 'Eukaryota genomes', y = 'Scaffold N50 (Mb)') +
  theme_bw()

human_genome_scfN50_plot = ggplot(subset(euk_genomes, Organism.Name == 'Homo sapiens')) + 
  geom_boxplot(aes(Assembly.Submission.Year, Assembly.Stats.Scaffold.N50/1000000)) +
  labs(title = 'Homo sapiens', y = 'Scaffold N50 (Mb)') +
  theme_bw()

(euk_genome_contigN50_plot | human_genome_contigN50_plot) /
  (euk_genome_scfN50_plot | human_genome_scfN50_plot)

# year ~ assembly size, dot plot log scale
ggplot(euk_genomes) + 
  geom_point(aes(Assembly.Stats.Total.Sequence.Length/1000000, Assembly.Stats.Contig.N50, 
                 color = Assembly.Submission.Year), size = 0.1) +
  scale_y_log10() + scale_x_log10() + 
  labs(title = 'Eukaryota genomes', y = 'Contig N50 (Mb)', x = 'Assembly size (Mb)') +
  theme_bw()

ggplot(subset(euk_genomes, Organism.Name == 'Homo sapiens')) + 
  geom_point(aes(Assembly.Stats.Total.Sequence.Length/1000000, Assembly.Stats.Contig.N50, 
                 color = Assembly.Submission.Year), size = 0.5) +
  labs(title = 'Homo sapiens', y = 'Contig N50 (Mb)', x = 'Assembly size (Mb)') +
  theme_bw()


euk_genome_size_plot_ori =  ggplot(euk_genomes, aes(Assembly.Submission.Year, Assembly.Stats.Total.Sequence.Length/1000000)) + 
  geom_boxplot() + 
  labs(title = 'Eukaryota genomes', y = 'Assembly size (Mb)') +
  theme_bw()

euk_genome_size_plot_log =  ggplot(euk_genomes, aes(Assembly.Submission.Year, Assembly.Stats.Total.Sequence.Length/1000000)) + 
  geom_boxplot() + 
  scale_y_log10() + 
  labs(title = 'Eukaryota genomes', y = 'Assembly size (Mb)') +
  theme_bw()

euk_genome_size_plot_ori / euk_genome_size_plot_log


# year ~ zebra fish contig N50
ggplot(subset(euk_genomes, Organism.Name == 'Danio rerio')) +
  geom_point(aes(Assembly.Submission.Year, Assembly.Stats.Contig.N50))

model_list = c('Homo sapiens', 'Saccharomyces cerevisiae', 'Arabidopsis thaliana', 'Mus musculus',
  'Drosophila melanogaster', 'Oryza sativa', 'Zea mays', 'Triticum aestivum', 'Hordeum vulgare',
  'Oryza sativa Indica Group', 'Gallus gallus', 'Zea mays subsp. mays', 'Ovis aries', 'Hordeum vulgare subsp. vulgare')

# year ~ assembly size
ggplot(euk_genomes, aes(Assembly.Submission.Year, Assembly.Stats.Total.Sequence.Length/1000000)) + 
  geom_point() + 
  labs(title = 'Eukaryota genomes', y = 'Assembly size (Mb)') +
  theme_bw()

# year ~ assembly size (non model species)
ggplot(subset(euk_genomes, !Organism.Name %in% model_list), aes(Assembly.Submission.Year, Assembly.Stats.Total.Sequence.Length/1000000)) + 
  geom_boxplot() + 
  labs(title = 'Eukaryota genomes', y = 'Assembly size (Mb)') +
  theme_bw()

# year ~ contig N50 (genome size > 2Gb)
ggplot(subset(euk_genomes, Assembly.Stats.Total.Sequence.Length > 2000000000), 
       aes(Assembly.Submission.Year, Assembly.Stats.Contig.N50/1000000)) + 
  geom_boxplot() + 
  labs(title = 'Eukaryota genomes > 2Gb', y = 'Contig N50 (Mb)') +
  theme_bw()

# year ~ assembly size (model species)
ggplot(subset(euk_genomes, Organism.Name %in% model_list), aes(Assembly.Submission.Year, Assembly.Stats.Contig.N50/1000000)) + 
  geom_boxplot() + facet_wrap(~Organism.Name) + scale_y_log10() +
  labs(title = 'Eukaryota genomes', y = 'Contig N50 (Mb)') +
  theme_bw()


### plot top species
# maximum N50 ~ year
target_list
euk_genomes_max = subset(euk_genomes, Organism.Name %in% target_list) %>% 
  group_by(Organism.Name, Assembly.Submission.Year, .add = TRUE) %>%
  summarise(contig.N50_mean = max(Assembly.Stats.Contig.N50, na.rm=T))
euk_genomes_max$Assembly.Submission.Year = as.numeric(euk_genomes_max$Assembly.Submission.Year)
euk_genomes_max_plot = ggplot(euk_genomes_max, aes(Assembly.Submission.Year, contig.N50_mean/1000000, shape = Organism.Name)) + 
  geom_point() + 
  geom_line() + #facet_wrap(~Organism.Name) + 
  scale_y_log10() +
  labs(title = 'Eukaryota genomes', y = 'Max contig N50 (Mb)') +
  theme_bw()

