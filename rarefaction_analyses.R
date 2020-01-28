# Analysis Biowide: rarefaction analyses
# Manuscript: Multitaxon inventory reveals highly consistent biodiversity responses to ecospace variation
# Author: Tobias Guldberg Frøslev
# Date: 27-05-2019

#setup
#### packages ####
library(iNEXT)
library(here)
library(ggplot2)
library(lemon)
#library(ggpmisc)
#library(ggpubr)
library(vegan)

#### load fungal data ####
fun_otu_tab <- readRDS(here::here("data","samples_focussed_otu_tab.RDS"))
fun_otu_tab <- fun_otu_tab[,1:130]

inext_fun <- iNEXT(fun_otu_tab, q=0, datatype="abundance", size=NULL, endpoint=NULL, knots=50, se=TRUE, conf=0.95, nboot=50)

saveRDS(inext_fun, here::here("data","inext_fun.RDS"))

rare_curve_fun <- ggiNEXT(inext_fun, type=1, se=TRUE, facet.var="site", color.var="order", grey=FALSE)  

formula <- y ~ x
rare_curve_fun2 <- rare_curve_fun + facet_rep_wrap(~site, ncol=8) +  theme_bw() + theme(panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                     text=element_text(size=15),
                     axis.text = element_text(colour="black"),
                     panel.border = element_blank(),
                     axis.line = element_line(colour = 'black', size = 0.25),
                     strip.background = element_rect(colour="white", fill="white")) + xlab("Sequencing depth / read count") + ylab("OTU richness") 

ggsave(here::here("plots","inext_fun_accumulation.pdf"), plot = rare_curve_fun2, device = "pdf", limitsize = T, width = 20, height = 30)

#### load eukaryote data ####
euk_otu_tab <- readRDS(here::here("data","euk_raw_derep.RDS"))

inext_euk <- iNEXT(euk_otu_tab, q=0, datatype="abundance", size=NULL, endpoint=NULL, knots=50, se=TRUE, conf=0.95, nboot=50)

saveRDS(inext_euk, here::here("data","inext_euk.RDS"))

rare_curve_euk <- ggiNEXT(inext_euk, type=1, se=TRUE, facet.var="site", color.var="order", grey=FALSE)

formula <- y ~ x
rare_curve_euk2 <- rare_curve_euk + facet_rep_wrap(~site, ncol=8) +  theme_bw() + theme(panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                     text=element_text(size=15),
                     axis.text = element_text(colour="black"),
                     panel.border = element_blank(),
                     axis.line = element_line(colour = 'black', size = 0.25),
                     strip.background = element_rect(colour="white", fill="white")) + xlab("Sequencing depth / read count") + ylab("OTU richness") 

ggsave(here::here("plots","inext_euk_accumulation.pdf"), plot = rare_curve_euk2, device = "pdf", limitsize = T, width = 20, height = 30)

#### load arthropod data ####
art_otu_tab <- readRDS(here::here("data","malaise_arthropod_combined.RDS"))
art_otu_tab <- art_otu_tab[,1:260]

inext_art <- iNEXT(art_otu_tab, q=0, datatype="abundance", size=NULL, endpoint=NULL, knots=50, se=TRUE, conf=0.95, nboot=50)

saveRDS(inext_art, here::here("data","inext_art.RDS"))

rare_curve_art <- ggiNEXT(inext_art, type=1, se=TRUE, facet.var="site", color.var="order", grey=FALSE)  

formula <- y ~ x
rare_curve_art2 <- rare_curve_art + facet_rep_wrap(~site, ncol=8) +  theme_bw() + theme(panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                     text=element_text(size=15),
                     axis.text = element_text(colour="black"),
                     panel.border = element_blank(),
                     axis.line = element_line(colour = 'black', size = 0.25),
                     strip.background = element_rect(colour="white", fill="white")) + xlab("Sequencing depth / read count") + ylab("OTU richness") 

ggsave(here::here("plots","inext_art_accumulation.pdf"), plot = rare_curve_art2, device = "pdf", limitsize = T, width = 20, height = 30)
