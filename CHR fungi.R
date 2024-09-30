packages <- c('ggcorrplot','ggthemes','dplyr', "ape", "ShortRead", "Biostrings",
              "phyloseq", "corncob","DESeq2", "microbiome", "DECIPHER", "phangorn",
              "tibble", "lme4", "lmerTest", "ggplot2", "vegan", "car", "rcompanion", 'microbiomeSeq',
              "emmeans", "RVAideMemoire",'gplots','plotly','tidyr','VennDiagram','venneuler')
sapply(packages, require, character.only = TRUE)              

# Set working directory
setwd("G:/My Drive/labs/NCSU/Projecto Shaneka 2/Fungi/")


# taxonomy ####
taxaCom_F <- readRDS("G:/My Drive/labs/NCSU/Projecto Shaneka 2/Fungi/tax_final240_240.rds")

###ASV table
asv.tableComF<- readRDS('G:/My Drive/labs/NCSU/Projecto Shaneka 2/Fungi/seqtab_final240_240.rds')
asv.tableComF<- otu_table(asv.tableComF, taxa_are_rows=FALSE)
rownames(asv.tableComF)

###Metadata
metac2 <- read.table("Metadata.txt", header = TRUE, row.names = 1)

#Now we can make the phyloseq object
ps.CompF <- phyloseq(asv.tableComF, tax_table(taxaCom_F), sample_data(metac2))
ps.CompF2 = subset_taxa(ps.CompF, Kingdom  == "k__Fungi")
taxa_names(ps.CompF2) <- paste0("ASV", seq(ntaxa(ps.CompF2)))
ps.CompF3 <- prune_samples(sample_sums(ps.CompF2)>=5, ps.CompF2)
ps.CompF4 <- prune_taxa(taxa_sums(ps.CompF3)>0, ps.CompF3)

ps.perc.F <- microbiome::transform(ps.CompF4, "compositional")
ps.H <- transform(ps.CompF4, "hellinger")

ps.perc.F@otu_table
ps.H@otu_table

ordered(sample_sums(ps.CompF4))

# 1. NMDS plot####
paleta_alive <- c("#C200FF",'#FF0000','#8B7500','#00008B',"#FFB919",'#FF7F50',"#00CC7A", 'black', 'grey')
#NRO shadeland Good

plot_ordination(ps.H, ordinate(ps.H, "PCoA", "bray"), color = "Location", shape = 'Species') + theme_few() +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  geom_point(size = 4) +
  scale_color_manual(values = paleta_alive)


#####PCoA for ASV-level data with Bray-Curtis
NMDS <- ordinate(ps.CompF.perc, "NMDS", "bray")
NMDS$points
x<-data.frame(NMDS$points[, 1:2])
x <- merge(x,ps.CompF.perc@sam_data, by=0)
rownames(x)<-x$Row.names
x$Row.names <- NULL


PCoAmean = aggregate(cbind(mean.x=MDS1, mean.y=MDS2)~Species*State,x,mean)
PCoAds = aggregate(cbind(ds.x=MDS1,ds.y=MDS2)~Species*State,x,sd)
PCoAmean2 <- cbind(PCoAmean, PCoAds[,3:4])

ggplot(x, aes(MDS1, MDS2)) + 
  geom_point(aes(color=Species, shape = State, size=3))

ggplot(data = PCoAmean2, aes(mean.x, mean.y)) +
  geom_vline(xintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_errorbarh(mapping = aes(xmin =mean.x - ds.x ,xmax = mean.x + ds.x),size = 0.1, alpha = 0.5) +
  geom_errorbar(mapping = aes(ymin =mean.y - ds.y ,ymax = mean.y + ds.y),size = 0.1, alpha = 0.5)  + 
  geom_point( size = 3, aes(col = Species, shape= State),stroke = 1) + 
  scale_color_manual(values = paleta_alive) +
  labs(y='NMDS2', x='NMDS1')+
  theme_few()
  
# 2. PERMANOVA
# Calculate bray curtis distance matrix
ps.perc_bray <- phyloseq::distance(ps.CompF.perc, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(ps.CompF.perc))
# Adonis test
PERMANOVA = adonis2(ps.perc_bray ~ Species*Location, data = sampledf)
PERMANOVA2 = data.frame(PERMANOVA)
variable  = c('variable','variable','variable','variable','variable')
PERMANOVA2$Significance <- "No Significant"
pval_thres <- 0.05
PERMANOVA2$Significance[which(PERMANOVA2$Pr..F. < pval_thres)] <- "Significant"
PERMANOVA2$Significance <- PERMANOVA2$Significance %>% factor

LC = data.frame(cbind(variable, row.names(PERMANOVA2), PERMANOVA2[,c(3,6)]))
colnames(LC) = c('variable','Effect', 'R2','Significance')
LC$Effect[which(LC$Effect == 'Species:State')] <- "Species:Location"
LC$Effect[which(LC$Effect == 'State')] <- "Location"


#Plot Variance
dim(LC)
colnames(LC)
dim(LC)
data_melt <- LC[-(5),]

#Reorder x axis: the order of the factor will be the same as in the data.csv file
data_melt$Effect <- as.character(data_melt$Effect)#Turn your 'treatment' column into a character vector
data_melt$Effect <- factor(data_melt$Effect, levels=unique(rev(data_melt$Effect)))#Then turn it back into a factor with the levels in the correct order

mypal = c('white','#dcda65','#528501','#1f4f58')

ggplot(data=data_melt, aes(x=variable, y=R2,fill=Effect)) +
  geom_bar(stat="identity",aes(color=Significance), size = 0.3,width = 1,linewidth=0.4) +# 
  scale_fill_manual(values = mypal) +
  theme_few() + #guides() + #color = FALSE, fill=FALSE
  labs(y="Variance explained (%)") +
  scale_color_manual(values = c('black',"red"),na.value =  "transparent",name = "Significance vs Water")+
  theme(legend.title=element_blank(), legend.margin=margin(c(0,0,0,0)),
        axis.text.x   = element_blank())
ggsave("Bacteria_variance.pdf", height=6, width=3, units='in')


# 3. Alfadiversity
# let's calcuolate Shannon diversity and Richness (Observed)
alpha.div <- estimate_richness(ps.CompF4, measures=c("Shannon", "Observed",'Chao', 'Simpson'))
even <- evenness(ps.1, 'pielou')
write.table(alpha.div, "Alfa fungi.txt")

plot_richness(ps.8, x = "Species", color = 'State',
              title = 'Alphadiversity', scales = "free_y", nrow = 1,
              measures = c("Shannon"), sortby = NULL)

Metadata <- ps.CompF4@sam_data

# it is easiest to append the alpha-diversity estimates to the metadata file for downstream analyses
Metadata$Shannon <- paste(alpha.div$Shannon)
Metadata$Observed <- paste(alpha.div$Observed)
Metadata$Simpson <- paste(alpha.div$Simpson)
Metadata$Chao1 <- paste(alpha.div$Chao1)
Metadata$Observed <- as.numeric(Metadata$Observed)
Metadata$Shannon <- as.numeric(Metadata$Shannon)
Metadata$Simpson <- as.numeric(Metadata$Simpson)
Metadata$Chao1 <- as.numeric(Metadata$Chao1)
Metadata$Location <- as.factor(Metadata$Location)
Metadata$Species <- as.factor(Metadata$Species)
Metadata$Replicates <- as.factor(Metadata$Replicates)
Metadata$Field <- as.factor(Metadata$Field)

#write.table(Metadata, "Alfadiversity fun.txt")


###Table ALphadiv
##Figure Alfa Diversidad Compartment

Metadata <- read.table("Alfadiversity fun.txt", header = TRUE, row.names = 1)
Metadata2 = Metadata[Metadata$Shannon > 0.5,]

ggplot(Metadata2, aes(Species, Shannon,  fill = State) ) +
  geom_boxplot() +
  theme_few()+
  labs(title ="Shannon index")+
  theme(legend.position="left")+
  ylim(c(0,NA))+
  #facet_grid(~State, scales = "free", space = "free")+
  scale_fill_manual(values = safe_colorblind_palette) +
  scale_color_manual(values = safe_colorblind_palette)


summary(Observed.model<-lmer(Shannon ~ State*Species  + (1|Replicates), data = Metadata))
anova(Observed.model) 
#Type III Analysis of Variance Table with Satterthwaite's method
#              Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
#State         2.6621  1.3311     2    24  10.481 0.0005349 ***
#Species       3.4789  3.4789     1    24  27.395 2.301e-05 ***
#State:Species 0.1397  0.0698     2    24   0.550 0.5840683    



# 4. Plotting Relative Abundance Bar Charts####
# phylum-level

ps.phyla <- tax_glom(ps.CompF4, "Genus")
ps.CompF.perc<- transform_sample_counts(ps.phyla, function(x) x / sum(x))
ps.CompF.perc@otu_table
ps.CompF.perc@tax_table
ordered(sample_sums(ps.CompF.perc))

write.table(ps.CompF4@otu_table, "otu table total fungi.txt")

# identify the 10 most abundant phylum
phyla = cbind(ps.CompF.perc@sam_data, ps.CompF.perc@otu_table)

melt.phylum <- melt(Phylum[,-c(2,4,7)])
head(melt.phylum)
safe_colorblind_palette <- c("#12AA99", "#882255", "#999933", "#6699CC", "#888888","#AB1499",
                             "#44AA09", "#999933", "#882255", "#661100", "#6699CC", "#886888")



melt.phylum$Treatment <- paste(melt.phylum$Species, melt.phylum$State, sep="_")
melt.phylum_mean <- melt.phylum%>%group_by(Species, variable)%>%
  summarise_all(mean)

ggplot(melt.phylum_mean, aes(x = Species, y = value, fill = variable)) + theme_few() +
  geom_bar(stat = "identity") + scale_fill_manual(values = safe_colorblind_palette) + labs(fill = "Phylum")+
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5))+ 
  ylab("Relative abundance")



melt.phylum_mean <- melt.phylum%>%group_by(State, sample_Species, Species)%>%
  summarise_all(mean)
melt.phylum_mean2 = data.frame(melt.phylum_mean)

ggplot(data=melt.phylum_mean2, aes(x=sample_Species, y=Abundance, fill=Species)) + 
  geom_bar(stat="identity")+ 
  ylab("Relative abundance (%)")+
  facet_grid(~State, scales = "free", space = "free")+
  scale_fill_manual(values = safe_colorblind_palette) +
  theme_few()

####Phylum analysis
ps.phyla
ncol(ps.phyla@otu_table)
Phylum = data.frame(cbind(ps.CompF.perc@sam_data,ps.CompF.perc@otu_table))
colnames(Phylum) = c("ID","Field","Species","Year_planted","Location","State","Replicates",
                                     "p__Basidiomycota","p__Ascomycota","p__Chytridiomycota","p__Mortierellomycota",
                                     "p__Glomeromycota")


#p__Basidiomycota
summary(aov(p__Basidiomycota ~ Species*Location   + (1|Replicates), data = Phylum))

#Df Sum Sq Mean Sq F value   Pr(>F)    
#Species        1  3.331   3.331  38.534 2.05e-06 ***
#State          2  0.003   0.001   0.015    0.985    
#Species:State  2  0.034   0.017   0.197    0.823    
#Residuals     24  2.075   0.086 

#p__Ascomycota
summary(aov(p__Ascomycota ~ Species*Location   + (1|Replicates), data = Phylum))
#Df Sum Sq Mean Sq F value Pr(>F)    
#Species        1  3.424   3.424  42.300  1e-06 ***
#State          2  0.026   0.013   0.161  0.852    
#Species:State  2  0.079   0.040   0.490  0.618    
#Residuals     24  1.943   0.081                   

#p__Chytridiomycota
summary(aov(p__Chytridiomycota ~ Species*Location   + (1|Replicates), data = Phylum))
#ns

#p__Mortierellomycota
summary(aov(p__Mortierellomycota ~ Species*Location   + (1|Replicates), data = Phylum))
#Df  Sum Sq  Mean Sq F value Pr(>F)  
#Species        1 0.00058 0.000577   0.131 0.7201  
#State          2 0.01437 0.007184   1.637 0.2156  
#Species:State  2 0.04080 0.020398   4.647 0.0197 *
#Residuals     24 0.10535 0.004390      

#p__Glomeromycota
summary(aov(p__Glomeromycota ~ Species*Location   + (1|Replicates), data = Phylum))
#ns

Phylum_mean <- Phylum%>%  summarise_all(mean)
Phylum_mean <- data.frame(Phylum_mean[,-c(2:7)])
#p__Basidiomycota p__Ascomycota p__Chytridiomycota p__Mortierellomycota
#0.4763619        0.4844196     0.009107981           0.03000277
#p__Glomeromycota 
#0.0001077006

#Alfadiversity and soil parameters corr
library(ggcorrplot)
AlfaTable = merge(Metadata, Table[,-c(1,3:6)], by='ID')

tablaCor<-cor(AlfaTable[,-(1:9)], method = 'spearman')
p.mat <-cor_pmat(AlfaTable[,-(1:9)], method = 'spearman')
ggcorrplot(tablaCor, p.mat = p.mat, hc.order = TRUE,
           type = 'lower', insig='blank')

#Phylla and soil parameters corr
library(ggcorrplot)
PhylumTable = merge(Phylum, Table[,-c(1,3:6)], by='ID')

tablaCor<-cor(PhylumTable[,-c(1:7,13)], method = 'spearman')
p.mat <-cor_pmat(PhylumTable[,-c(1:7,13)], method = 'spearman')
ggcorrplot(tablaCor, p.mat = p.mat, hc.order = TRUE,
           type = 'lower', insig='blank')



#Species-level

ps.phyla <- tax_glom(ps.CompF4, "Species")
ps.CompF.perc<- transform_sample_counts(ps.phyla, function(x) x / sum(x))

species_tax = data.frame(ps.CompF.perc@tax_table)
species_tax$genus_species <- paste(species_tax$Genus, species_tax$Species, sep="_")
Species_ID = species_tax$genus_species 

Species = cbind(ps.CompF.perc@sam_data[,c(3,6)], ps.CompF.perc@otu_table)
colnames(Species) = c('Species','State',"Russula farinipes","Tuber linsdalei","Cortinarius saniosus","Trichoderma lixii",
                      "Achroceratosphaeria potamia", "Archaeorhizomyces sborealis","Atrocalyx nordicus","Russula pulverulenta",
                      "Russula samoenolens","Russula recondita","Tomentella tenuirhizomorpha","Archaeorhizomyces finlayi",
                      "Woswasia atropurpurea","Talaromyces rufus","Metarhizium brunneum","Densocarpa shanorii",        
                      "Trichoderma pubescens","Phallus rugulosus","Podila minutissima","Coniochaeta verticillata" )

melt.Species <- melt(Species)
melt.Species2 <- na.omit (melt.Species, "NanN") 

log2foldchange <- c()
log2foldchange_Complete<- c()

for(site in melt.Species2$State%>% unique){
  melted_sub2 <- melt.Species2 %>% subset(State ==  site) %>% droplevels
  for(specie in melted_sub2$variable  %>% unique){
    melted_sub <- melted_sub2 %>% subset(variable ==  specie) %>% droplevels #'ASV522'
    
    BW = data.frame(t(melted_sub[melted_sub$Species=='BW',][,4]))
    NR = data.frame(t(melted_sub[melted_sub$Species=='NRO',][,4]))
    
    BW_mean = apply(BW, 1, mean) 
    NR_mean = apply(NR, 1, mean) 
    
    BWmean_NRmean <- log2(BW_mean+1) - log2(NR_mean+1) 
    
    BWmean_NRmean_statistic = t.test(BW,NR)
    pvalue = BWmean_NRmean_statistic$p.value#Water_mean+1
    
    result = cbind(BWmean_NRmean,pvalue)
    result2 = cbind(result,site)
    row.names(result2)= specie
    log2foldchange <- rbind(log2foldchange,result2)
  }
  log2foldchange_Complete = data.frame(cbind(row.names(log2foldchange),log2foldchange))
  colnames(log2foldchange_Complete)=c('Species','diff','pvalue','Location')
}

log2foldchange_Complete$Significance <- "No Significant"
pval_thres <- 0.05
log2foldchange_Complete$Significance[which(log2foldchange_Complete$pvalue < pval_thres)] <- "q < 0.05"
log2foldchange_Complete$Significance <- log2foldchange_Complete$Significance %>% factor

log2foldchange_Complete$diff = as.numeric(log2foldchange_Complete$diff)

ggplot(data = log2foldchange_Complete, aes(Location,Species)) +
  geom_raster(aes(fill = diff)) +
  geom_tile(aes(color = Significance),fill = '#00000000', size = 0.3,width = 0.9,height = 0.95) + #
  #facet_grid(~Compartment, space = "free",scales = "free") +
  theme_few()+
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_40_95_c42", #kovesi.diverging_bwr_55_98_c37
                         limits = c(-1,1),na.value = "#D9D9D9",name = "Fold Change") +
  scale_color_manual(values = c('black',"red"),na.value =  "transparent",name = "Significance") + #Significance Genotype vs Col-0
  theme(axis.text.x = element_text(angle = -45, hjust=-0.05, size = 4),
        axis.title = element_blank())
