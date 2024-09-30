library("FactoMineR")
library("factoextra")
library(ade4)
library(MASS)
library(ellipse)
library(ggplot2)
library(agricolae)
library(vegan)
library(ggthemes)
library(RVAideMemoire)

# Set working directory
setwd("C:/Users/juanp/Google Drive/labs/NCSU/Projecto Shaneka 2/")

Table <- read.table("Soil properties.txt", header = TRUE)
head(Table)

Table2 = Table[Table$ID!='B212',]

melt.Table <- Table3%>%group_by(Plant_state, State, Species)%>%
  summarise_all(mean)

summary(lm(Height~Dbh, data=Table3))

Table3$Plant_state = factor(Table3$Plant_state, c('Good','Bad'))
#Height
Table3 = Table2[Table2$Height < 100,]
melt.Table[,c(1:4,36,35)]

aovHeight <- aov(Height~Plant_state*State*Species, data=Table3)
summary(aovHeight)
#Plant_state                1  624.6   624.6 279.459 < 2e-16 ***
#State                      2  469.2   234.6 104.948 < 2e-16 ***
#Species                    1    0.5     0.5   0.240 0.62585    
#Plant_state:State          2   30.7    15.3   6.866 0.00175 ** 
#Plant_state:Species        1   15.1    15.1   6.734 0.01120 *  
#State:Species              2   30.9    15.4   6.910 0.00169 ** 
#Plant_state:State:Species  2   12.4     6.2   2.771 0.06847 .  
#Residuals                 82  183.3     2.2    

lsdHeight <- LSD.test(aovHeight, c('Plant_state','Species'), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdHeight


safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499") 


ggplot(data=Table3, aes(x=Species, y=Height, fill=Plant_state)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~State,space = "fixed",scales = "free_y") +
  ggtitle("Height")+ ylab("meters")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Plant Species")+ theme_bw()+ylim(0,NA)

#Dbh
aovDbh <- aov(Dbh~Plant_state*State*Species, data=Table)
summary(aovDbh)
#Plant_state                  1 2338.2  2338.2 502.936 < 2e-16 ***
#State                        2 1020.1   510.1 109.710 < 2e-16 ***
#Species                      1   32.9    32.9   7.075 0.00938    
#Plant_state:State            2   14.1     7.1   1.521  0.22458    
#Plant_state:Species          1   81.6    81.6  17.560 6.91e-05 ***
#State:Species                2  219.0   109.5  23.549 7.93e-09 ***
#Plant_state:State:Species    2   30.8    15.4   3.313  0.04126 *  
#Residuals                    3  385.9     4.6       

lsdDbh <- LSD.test(aovDbh, c('Plant_state','Species'), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdDbh


safe_colorblind_palette <- c("#DDCC77","#CC6677", "#88CCEE",  "#117733", "#332288", "#AA4499") 

ggplot(data=Table3, aes(x=Species, y=Dbh, fill=Plant_state)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~State,space = "fixed",scales = "free_y") +
  ggtitle("Dbh")+ ylab("meters")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Plant Species")+ theme_bw()+ylim(0,NA)

#Mg
aovMg <- aov(Mg~Plant_state*State*Species, data=Table2)
summary(aovMg)

#Plant_state                1   3.72   3.720   4.919   0.0293 *  
#State                      2   5.47   2.736   3.617   0.0312 *  
#Species                    1  26.48  26.480  35.011 7.30e-08 ***
#Plant_state:State          2   0.34   0.170   0.225   0.7987    
#Plant_state:Species        1   0.80   0.800   1.058   0.3068    
#State:Species              2  35.01  17.504  23.144 1.07e-08 ***
#Plant_state:State:Species  2   0.88   0.441   0.583   0.5606    
#Residuals                 82  62.02   0.756   

lsdMg <- LSD.test(aovMg, c('Plant_state','State','Species'), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdMg

summary(lm(Height~Mg, data=Table3))
ggplot(data=Table3, aes(x=Species, y=Height, fill=Plant_state)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~State,space = "fixed",scales = "free_y") +
  ggtitle("Height")+ ylab("meters")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Plant Species")+ theme_bw()+ylim(0,NA)

Table3$Mg <- as.numeric(Table3$Mg)

summary(lm(Height~Mg, data=Table3))

ggplot(Table3, aes(y=Height, x=Mg)) + 
  geom_point(aes(col=Plant_state)) + 
  geom_smooth(method="loess", se=T,aes(col=Plant_state))


ggplot(data=Table2, aes(x=Species, y=Mg, fill=Plant_state)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~State,space = "fixed",scales = "free_y") +
  ggtitle("Mg")+ ylab("meters")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Plant Species")+ theme_bw()+ylim(0,NA)


#Ca
aovCa <- aov(Ca~Plant_state*State*Species, data=Table[Table$DNA!='212',])
summary(aovCa)

#Plant_state                1   15.2   15.20   1.953   0.1660    
#State                      2   42.7   21.36   2.746   0.0701 .  
#Species                    1  149.4  149.37  19.196 3.46e-05 ***
#Plant_state:State          2    3.1    1.55   0.200   0.8193    
#Plant_state:Species        1   16.9   16.87   2.168   0.1447    
#State:Species              2   65.7   32.84   4.220   0.0180 *  
#Plant_state:State:Species  2    8.0    4.01   0.515   0.5995    
#Residuals                 82  638.1    7.78                     

lsdCa <- LSD.test(aovCa, c('Plant_state','Species'), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdCa


Table[Table$DNA!='212',]$Ca <- as.numeric(Table[Table$DNA!='212',]$Ca)

ggplot(data=Table2, aes(x=Species, y=Ca, fill=Plant_state)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~State,space = "fixed",scales = "free_y") +
  ggtitle("Ca")+ ylab("meters")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Plant Species")+ theme_bw()+ylim(0,NA)


#Na
aovNa <- aov(Na~Plant_state*State*Species, data=Table2)
summary(aovNa)

#Plant_state                1 0.00011 0.000106   0.131  0.718  
#State                      2 0.00686 0.003429   4.218  0.018 *
#Species                    1 0.00182 0.001823   2.242  0.138  
#Plant_state:State          2 0.00105 0.000524   0.645  0.528  
#Plant_state:Species        1 0.00046 0.000460   0.565  0.454  
#State:Species              2 0.00335 0.001674   2.059  0.134  
#Plant_state:State:Species  2 0.00107 0.000536   0.660  0.520  
#Residuals                 82 0.06667 0.000813 


lsdNa <- LSD.test(aovNa, c('Plant_state','Species'), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdNa



Table2$Na <- as.numeric(Table2$Na)

ggplot(data=Table2, aes(x=Species, y=Na, fill=Plant_state)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~State,space = "fixed",scales = "free_y") +
  ggtitle("Na")+ ylab("meters")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Plant Species")+ theme_bw()+ylim(0,NA)

#K
aovK <- aov(K~Plant_state*State*Species, data=Table2)
summary(aovK)

#Plant_state                1 0.0086  0.0086   0.366 0.547071    
#State                      2 0.6971  0.3485  14.789 3.28e-06 ***
#Species                    1 0.0360  0.0360   1.529 0.219763    
#Plant_state:State          2 0.0798  0.0399   1.692 0.190436    
#Plant_state:Species        1 0.0008  0.0008   0.032 0.858512    
#State:Species              2 0.3783  0.1892   8.026 0.000656 ***
#Plant_state:State:Species  2 0.0392  0.0196   0.833 0.438587    
#Residuals                 82 1.9326  0.0236 

lsdK <- LSD.test(aovK, c('Plant_state','Species'), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdK



Table2$K <- as.numeric(Table2$K)

ggplot(data=Table2, aes(x=Species, y=K, fill=Plant_state)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~State,space = "fixed",scales = "free_y") +
  ggtitle("K")+ ylab("meters")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Plant Species")+ theme_bw()+ylim(0,NA)

#SumNH4
aovSumNH4 <- aov(SumNH4~Plant_state*State*Species, data=Table2)
summary(aovSumNH4)

#Plant_state                1   32.8   32.77   2.455  0.12097    
#State                      2   60.5   30.25   2.267  0.11009    
#Species                    1  309.8  309.77  23.211 6.56e-06 ***
#Plant_state:State          2    6.9    3.43   0.257  0.77369    
#Plant_state:Species        1   25.0   24.96   1.870  0.17522    
#State:Species              2  180.8   90.41   6.775  0.00189 ** 
#Plant_state:State:Species  2   14.5    7.26   0.544  0.58263    
#Residuals                 82 1094.3   13.35 


aovSumNH4 <- aov(SumNH4~Plant_state, data=Table2[Table2$Species=='NRO',])
summary(aovSumNH4)

lsdSumNH4 <- LSD.test(aovSumNH4, c('Plant_state'), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdSumNH4

##NRO hay algo ahi con el N.

safe_colorblind_palette <- c("#CC6677","#DDCC77", "#88CCEE",  "#117733", "#332288", "#AA4499") 

Table2$SumNH4 <- as.numeric(Table2$SumNH4)

ggplot(data=Table2, aes(x=Species, y=SumNH4, fill=Plant_state)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~State,space = "fixed",scales = "free_y") +
  ggtitle("SumNH4")+ ylab("meters")+ 
  theme(plot.title = element_text(color="blacSumNH4", size=17, face="bold.italic"))+ 
  xlab("Plant Species")+ theme_bw()+ylim(0,NA)



#CEC
aovCEC <- aov(CEC~Plant_state*State*Species, data=Table2)
summary(aovCEC)

#                          Df Sum Sq Mean Sq F value  Pr(>F)   
#Plant_state                1    6.9    6.86   0.667 0.41634   
#State                      2   17.5    8.74   0.850 0.43125   
#Species                    1   92.5   92.47   8.991 0.00359 **
#Plant_state:State          2    0.8    0.41   0.040 0.96104   
#Plant_state:Species        1   11.7   11.74   1.141 0.28848   
#State:Species              2   63.5   31.73   3.085 0.05106 . 
#Plant_state:State:Species  2    2.5    1.26   0.122 0.88508   
#Residuals                 82  843.3   10.28   


lsdCEC <- LSD.test(aovCEC, c('Plant_state','Species'), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdCEC


safe_colorblind_palette <- c("#CC6677","#DDCC77", "#88CCEE",  "#117733", "#332288", "#AA4499") 

Table2$CEC <- as.numeric(Table2$CEC)

ggplot(data=Table2, aes(x=Species, y=CEC, fill=Plant_state)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~State,space = "fixed",scales = "free_y") +
  ggtitle("CEC")+ ylab("meters")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Plant Species")+ theme_bw()+ylim(0,NA)

#Base_Sat
aovBase_Sat <- aov(Base_Sat~Plant_state*State*Species, data=Table2)
summary(aovBase_Sat)

#                         Df Sum Sq Mean Sq F value   Pr(>F)    
#Plant_state                1   1995    1995   9.306  0.00308 ** 
#State                      2   1471     735   3.432  0.03705 *  
#Species                    1   7510    7510  35.041 7.22e-08 ***
#Plant_state:State          2   1492     746   3.482  0.03537 *  
#Plant_state:Species        1    802     802   3.742  0.05652 .  
#State:Species              2   2539    1270   5.924  0.00395 ** 
#Plant_state:State:Species  2    682     341   1.591  0.20998    
#Residuals                 82  17574     214    

lsdBase_Sat <- LSD.test(aovBase_Sat, c('Plant_state','Species'), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdBase_Sat


Table2$Base_Sat <- as.numeric(Table2$Base_Sat)

ggplot(data=Table2, aes(x=Species, y=Base_Sat, fill=Plant_state)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~State,space = "fixed",scales = "free_y") +
  ggtitle("Base_Sat")+ ylab("meters")+ 
  theme(plot.title = element_text(color="blacBase_Sat", size=17, face="bold.italic"))+ 
  xlab("Plant Species")+ theme_bw()+ylim(0,NA)

Table3$Base_Sat <- as.numeric(Table3$Base_Sat)

summary(lm(Height~Base_Sat, data=Table3))

ggplot(Table3, aes(y=Height, x=Base_Sat)) + 
  geom_point() + 
  geom_smooth(method="loess", se=T) #,aes(col=Plant_state)

#Corg
aovCorg <- aov(Corg~Plant_state*State*Species, data=Table2)
summary(aovCorg)

#                          Df Sum Sq Mean Sq F value   Pr(>F)    
#Plant_state                1  0.013   0.013   0.066   0.7980    
#State                      2  0.008   0.004   0.022   0.9786    
#Species                    1  4.352   4.352  22.291 9.54e-06 ***
#Plant_state:State          2  0.055   0.028   0.141   0.8685    
#Plant_state:Species        1  0.001   0.001   0.004   0.9497    
#State:Species              2  0.941   0.470   2.410   0.0962 .  
#Plant_state:State:Species  2  0.136   0.068   0.347   0.7075    
#Residuals                 82 16.008   0.195                  

lsdCorg <- LSD.test(aovCorg, c('Plant_state','Species'), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdCorg


Table2$Corg <- as.numeric(Table2$Corg)

ggplot(data=Table2, aes(x=Species, y=Corg, fill=Plant_state)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~State,space = "fixed",scales = "free_y") +
  ggtitle("Corg")+ ylab("meters")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Plant Species")+ theme_bw()+ylim(0,NA)

#pH_CaCl2
aovpH_CaCl2 <- aov(pH_CaCl2~Plant_state*State*Species, data=Table2)
summary(aovpH_CaCl2)


#                          Df Sum Sq Mean Sq F value   Pr(>F)    
#Plant_state                1  1.716   1.716   5.596  0.02036 *  
#State                      2  1.710   0.855   2.789  0.06731 .  
#Species                    1 17.298  17.298  56.418 6.44e-11 ***
#Plant_state:State          2  1.377   0.688   2.245  0.11239    
#Plant_state:Species        1  0.203   0.203   0.662  0.41822    
#State:Species              2  4.520   2.260   7.370  0.00114 ** 
#Plant_state:State:Species  2  0.683   0.342   1.114  0.33304    
#Residuals                 82 25.141   0.307  



lsdpH <- LSD.test(aovpH_CaCl2, c('Plant_state','Species','State'), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdpH

Table2$pH_CaCl2 <- as.numeric(Table2$pH_CaCl2)

ggplot(data=Table2, aes(x=Species, y=pH_CaCl2, fill=Plant_state)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~State,space = "fixed",scales = "free_y") +
  ggtitle("pH CaCl2")+ ylab("meters")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Plant Species")+ theme_bw()+ylim(0,NA)

#pH_H2O
aovpH_H2O <- aov(pH_H2O~Plant_state*State*Species, data=Table2)
summary(aovpH_H2O)


#                          Df Sum Sq Mean Sq F value   Pr(>F)    
#Plant_state                1  1.383   1.383   5.267 0.024288 *  
#State                      2  0.885   0.443   1.686 0.191547    
#Species                    1 16.145  16.145  61.507 1.43e-11 ***
#Plant_state:State          2  1.063   0.531   2.024 0.138662    
#Plant_state:Species        1  0.210   0.210   0.801 0.373539    
#State:Species              2  4.495   2.248   8.563 0.000419 ***
#Plant_state:State:Species  2  0.652   0.326   1.242 0.294228    
#Residuals                 82 21.524   0.262 

lsdpH <- LSD.test(aovpH_H2O, c('Plant_state','Species','State'), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdpH

Table3$pH_H2O <- as.numeric(Table3$pH_H2O)

ggplot(data=Table2, aes(x=Species, y=pH_H2O, fill=Plant_state)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~State,space = "fixed",scales = "free_y") +
  ggtitle("pH H2O")+ ylab("meters")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Plant Species")+ theme_bw()+ylim(0,NA)

summary(lm(Height~pH_H2O, data=Table3[Table3$Species == "NRO",]))
summary(lm(Height~pH_H2O, data=Table3[Table3$Species == "BW",]))


ggplot(Table3, aes(y=Height, x=pH_H2O)) + 
  geom_point(aes(col=Species, shape=Plant_state)) + 
  facet_grid(.~State,space = "fixed",scales = "free") +
  geom_smooth(method="loess", se=T,aes(col=Species))+
  theme_bw()

ggplot(Table3[Table3$Species == "BW",], aes(y=Height, x=pH_H2O)) + 
  geom_point(aes(col=State, shape=Species, size=Height)) + 
  geom_smooth(method="loess", se=T,aes(col=State))+
  theme_bw()


#PBray
aovPBray <- aov(PBray~Plant_state*State*Species, data=Table2)
summary(aovPBray)


#                          Df Sum Sq Mean Sq F value   Pr(>F)    
#Plant_state                1    472   472.3   3.203 0.0772 .
#State                      2     23    11.4   0.078 0.9254  
#Species                    1    217   216.8   1.470 0.2288  
#Plant_state:State          2    689   344.6   2.337 0.1030  
#Plant_state:Species        1     49    49.0   0.332 0.5660  
#State:Species              2   1205   602.6   4.086 0.0203 *
#Plant_state:State:Species  2    711   355.5   2.410 0.0961 .
#Residuals                 82  12092   147.5


lsdpH <- LSD.test(aovPBray, c('Plant_state','Species','State'), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdpH

Table2$PBray <- as.numeric(Table2$PBray)

ggplot(data=Table2, aes(x=Species, y=PBray, fill=Plant_state)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~State,space = "fixed",scales = "free_y") +
  ggtitle("pH H2O")+ ylab("meters")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Plant Species")+ theme_bw()+ylim(0,NA)


#Sand
aovSand <- aov(Sand~Plant_state*State*Species, data=Table2)
summary(aovSand)


#                          Df Sum Sq Mean Sq F value   Pr(>F)    
#Plant_state                1    472   472.3   3.203 0.0772 .
#State                      2     23    11.4   0.078 0.9254  
#Species                    1    217   216.8   1.470 0.2288  
#Plant_state:State          2    689   344.6   2.337 0.1030  
#Plant_state:Species        1     49    49.0   0.332 0.5660  
#State:Species              2   1205   602.6   4.086 0.0203 *
#Plant_state:State:Species  2    711   355.5   2.410 0.0961 .
#Residuals                 82  12092   147.5


lsdpH <- LSD.test(aovSand, c('Plant_state','Species','State'), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdpH

Table2$Sand <- as.numeric(Table2$Sand)

ggplot(data=Table2, aes(x=Species, y=Sand, fill=Plant_state)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~State,space = "fixed",scales = "free_y") +
  ggtitle("pH H2O")+ ylab("meters")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Plant Species")+ theme_bw()+ylim(0,NA)

Table3$Sand <- as.numeric(Table3$Sand)

summary(lm(Height~Sand, data=Table3[Table3$Species == "NRO",]))
summary(lm(Height~Sand, data=Table3[Table3$Species == "BW",]))


ggplot(Table3, aes(y=Height, x=Sand)) + 
  geom_point(aes(col=State, shape=Species)) + 
  geom_smooth(method="loess", se=T,aes(col=Plant_state))+
  theme_bw()

ggplot(Table3, aes(y=Height, x=Sand)) + 
  geom_point() + 
  geom_smooth(method="loess", se=T) #

#Clay
aovClay <- aov(Clay~Plant_state*State*Species, data=Table2)
summary(aovClay)


#                          Df Sum Sq Mean Sq F value   Pr(>F)    
#Plant_state                1   25.3   25.33   1.159  0.285
#State                      2   99.7   49.85   2.280  0.109
#Species                    1    5.5    5.48   0.251  0.618
#Plant_state:State          2   10.0    5.01   0.229  0.796
#Plant_state:Species        1   12.8   12.79   0.585  0.447
#State:Species              2   45.8   22.92   1.048  0.355
#Plant_state:State:Species  2   11.2    5.62   0.257  0.774
#Residuals                 82 1793.1   21.87

lsdpH <- LSD.test(aovSand, c('Plant_state','Species','State'), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdpH

Table2$Sand <- as.numeric(Table2$Sand)

ggplot(data=Table2, aes(x=Species, y=Clay, fill=Plant_state)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~State,space = "fixed",scales = "free_y") +
  ggtitle("Clay")+ ylab("meters")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Plant Species")+ theme_bw()+ylim(0,NA)

#Silt
aovSilt <- aov(Silt~Plant_state*State*Species, data=Table2)
summary(aovSilt)


#                          Df Sum Sq Mean Sq F value   Pr(>F)    
#Plant_state                1     31      31   0.312 0.5781    
#State                      2  35431   17716 179.368 <2e-16 ***
#Species                    1     83      83   0.838 0.3626    
#Plant_state:State          2    790     395   3.997 0.0221 *  
#Plant_state:Species        1    190     190   1.927 0.1689    
#State:Species              2    531     265   2.686 0.0742 .  
#Plant_state:State:Species  2    180      90   0.913 0.4054    
#Residuals                 82   8099      99   


lsdpH <- LSD.test(aovSand, c('Plant_state','Species','State'), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdpH

Table2$Sand <- as.numeric(Table2$Sand)

ggplot(data=Table2, aes(x=Species, y=Silt, fill=Plant_state)) +
  geom_boxplot()+
  scale_fill_manual(values=safe_colorblind_palette)+
  facet_grid(.~State,space = "fixed",scales = "free_y") +
  ggtitle("Silt")+ ylab("meters")+ 
  theme(plot.title = element_text(color="black", size=17, face="bold.italic"))+ 
  xlab("Plant Species")+ theme_bw()+ylim(0,NA)

##Corr nutrients
TableNRO = data.frame(Table3[Table3$Species=='NRO',])
TableBW = Table2[Table2$Species!='NRO',]
Table2 = as.numeric(Table2[,19:32]$Ca)


Table2$Mg = as.numeric(Table2$Mg)
Table2$Ca = as.numeric(Table2$Ca)
Table2$K = as.numeric(Table2$K)
Table2$Na = as.numeric(Table2$Na)
Table2$SumNH4 = as.numeric(Table2$SumNH4)
Table2$Corg = as.numeric(Table2$Corg)
Table2$CEC = as.numeric(Table2$CEC)
Table2$Base_Sat = as.numeric(Table2$Base_Sat)
Table2$pH_CaCl2 = as.numeric(Table2$pH_CaCl2)
Table2$pH_H2O = as.numeric(Table2$pH_H2O)
Table2$PBray = as.numeric(Table2$PBray)
Table2$Clay = as.numeric(Table2$Clay)
Table2$Silt = as.numeric(Table2$Silt)
Table2$Sand = as.numeric(Table2$Sand)
Table2$Height = as.numeric(Table2$Height)
Table2$Dbh = as.numeric(Table2$Dbh)
ncol(Table2)
is.numeric(Table2$Dbh)

Corr_nut<-cor(Table3[,c(15:18,35,36)], method = 'pearson')
p.mat_nut <-cor_pmat(Table3[,c(15:18,35,36)], method = 'pearson')
ncol(Nutrients)
melted_Corr_nut <-Corr_nut %>% melt
melted_p.mat_nut <-p.mat_nut %>% melt
melted_nut = cbind(melted_Corr_nut, melted_p.mat_nut[,3])
colnames(melted_nut) = c('Var1','Var2','value', 'pvalue')
nrow(melted_nut)

melted_nut$Significance <- "NoSignificant"
pval_thres <- 0.05
melted_nut$Significance[which(melted_nut$pvalue < pval_thres)] <- "Significant"
melted_nut$Significance <- melted_nut$Significance %>% factor

dend_nut <- as.dendrogram(hclust(dist(Corr_nut)))
dend_nut_data <- dendro_data(dend_nut)
dend_nut_data_order = as.character(dend_nut_data$labels[,3])

nrow(melted_nut)

melted_nut$Var2 = factor(melted_nut$Var2, c('CEC','Ca','SumNH4','Mg','Base_Sat','pH_CaCl2','pH_H2O','Corg','Clay','K','PBray','Na','Silt','Sand'))
melted_nut$Var1 = factor(melted_nut$Var1, c('CEC','Ca','SumNH4','Mg','Base_Sat','pH_CaCl2','pH_H2O','Corg','Clay','K','PBray','Na','Silt','Sand'))


ggplot(data = melted_nut[c(225:238,241:254),], aes(Var1,Var2)) + #
  geom_raster(aes(fill = value))+
  theme_few() + 
  geom_tile(aes(color = Significance),fill = '#00000000', size = 0.2,width = 0.9,height = 0.95) + #
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-1,1),na.value = "#D9D9D9",name = "Corr") +
  scale_color_manual(values = c('grey',"black"),na.value =  "transparent",name = "Significative") +
  labs(y='Microbial paramenters',x='Chemical parameters')+
  theme(axis.text.x = element_text(angle = -45, hjust=-0.05),
        axis.title = element_blank())

###Random forest 
library(randomForest)

# Split into Train and Validation sets
# Training Set : Validation Set = 70 : 30 (random)
set.seed(100)
train <- sample(nrow(TableBW), 0.7*nrow(TableBW), replace = FALSE)
TrainSet1 <- TableBW[train,]
ValidSet1 <- TableBW[-train,]
summary(TrainSet1)
summary(ValidSet1)

# Create a Random Forest model with default parameters
model <- randomForest(y = TableBW$Height, x = TableBW[,c('CEC','Ca','SumNH4','Mg','Base_Sat','pH_CaCl2','pH_H2O','Corg','Clay','K','PBray','Na','Silt','Sand')], 
                      data = TrainSet1, importance = TRUE, ntree = 500, mtry = 8, proximity = TRUE)
model
plot(model)

# Predicting on train set
predTrain1 <- predict(model, TrainSet1, type = "class")
# Checking classification accuracy
table(predTrain1, TrainSet1$Height) 

# Predicting on Validation set
predValid1 <- predict(model, ValidSet1, type = "class")
# Checking classification accuracy
mean(predValid1 == ValidSet1$Height)
table(predValid1,ValidSet1$Height)

# To check important variables
BW_importance  = importance(model) #         
varImpPlot(model)

BW_importance = data.frame(cbind(rownames(BW_importance),BW_importance)) 
colnames(BW_importance) = c('V1','BW_IncMSE', 'BW_IncNodePurity')
BW_importance$BW_IncMSE = as.numeric(BW_importance$BW_IncMSE)
BW_importance$V1 = factor(BW_importance$V1, c('CEC','PBray','Na','Mg','Corg','K','Clay','SumNH4','Ca','pH_H2O','Base_Sat','pH_CaCl2','Silt','Sand'))


ggplot(data=BW_importance, aes(x=V1, y=BW_IncMSE, fill=BW_IncMSE)) +
  geom_bar(stat="identity")+
  theme_few()+ labs(y= '%IncMSE',x= 'Soil Properties')+
  coord_flip()+ theme(legend.position = 'none')


B_importance = data.frame(cbind(rownames(B_importance),B_importance)) 
colnames(B_importance) = c('V1','B_IncMSE', 'B_IncNodePurity')
B_importance$B_IncMSE = as.numeric(B_importance$B_IncMSE)

B_importance$V1<-factor(B_importance$V1,rev(c("CCL2_diameter","Epidermis_diameter","cell_layer_n","Perycicle_diameter",'Max_diameter',"Arenquima_area_prom","Still_diameter","Endodermis_diameter","Exodermis_diameter", "CCL1_diameter","CCL3_diameter","Epidermis_2_diameter")))
B_importance$B_IncMSE <- factor(B_importance$B_IncMSE, levels=unique(rev(B_importance$B_IncMSE)))#Then turn it back into a factor with the levels in the correct order

ggplot(data=B_importance, aes(x=V1, y=B_IncMSE, fill=B_IncMSE)) +
  geom_bar(stat="identity") +
  #scale_fill_manual(values = mypal) +
  theme_few()+ labs(y='Root traits',
                    x='%IncMSE')+
  coord_flip()+ theme(legend.position = 'none')


all_importance = data.frame(cbind(rownames(all_importance),all_importance)) 
colnames(all_importance) = c('V1','all_IncMSE', 'all_IncNodePurity')
all_importance$all_IncMSE = as.numeric(all_importance$all_IncMSE)
all_importance$V1<-factor(all_importance$V1,rev(c("Epidermis_diameter","Perycicle_diameter","Endodermis_diameter","Still_diameter","CCL2_diameter",'Max_diameter',"Arenquima_area_prom","cell_layer_n","CCL1_diameter","Exodermis_diameter", "CCL3_diameter","Epidermis_2_diameter")))

ggplot(data=all_importance, aes(x=V1, y=all_IncMSE, fill=all_IncMSE)) +
  geom_bar(stat="identity")+
  theme_few()+ labs(y= 'Root traits',x= '%IncMSE')+
  coord_flip()+ theme(legend.position = 'none')




