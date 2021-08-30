#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
#  Alpha Diversity
#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

library(dplyr) # Data Manipulation
library(plyr) # Data Manipulation
library(lsmeans) # pairwise stats
library(ggplot2) # plotting
library(vegan) # Diversity Stats
library(car) # Model Assumption Checks
library(ggpubr) # Graphing Model Assumption Checks

Richness<-read.csv("ARMS_Data_Working.csv")
metadata<-read.csv("MetaData_ARMS.csv")

Richness2<-merge(metadata, Richness, by = "Sample")
Richness2$Temperature<-as.factor(Richness2$Temperature)
Richness2$CO2<-as.factor(Richness2$CO2)
Richness2$HeaderTank<-as.factor(Richness2$HeaderTank)
Richness2$Richness<-specnumber(Richness2[,8:ncol(Richness2)])
Richness2 <- Richness2 %>% dplyr::select(Richness, everything())

        #><><><><><><><><><><><><><><><>
        #  Table S5- Summary Stats
        #><><><><><><><><><><><><><><><>
            SummaryRichness<-Richness2 %>%
            group_by(Treatment) %>%
            get_summary_stats(Richness, type = "mean_sd")
        
#><><><><><><><><><><><><><><>
#  ANOVA Richness 
#><><><><><><><><><><><><><><>
            
# ANOVA Richness ~ Conditions
model  <- aov(Richness~CO2*Temperature/HeaderTank, data=Richness2)
              #<><><><><><><><><><><><><>
              # Table S6 - summary output 
              #<><><><><><><><><><><><><>
summary(model) 
              #<><><><><><><><><><><><><>
              # Table S7 - lsmean output 
              #<><><><><><><><><><><><><>
lsmeans(model, pairwise ~ CO2*Temperature) 

# Checking residuals for assumptions of normality and homogeneity
ggqqplot(residuals(model))
plot(model,1)
# Computing Shapiro-Wilk test of normality
shapiro.test(residuals(model))
# Checking for homogeneity of variance
bartlett.test(Richness ~ interaction(CO2,Temperature),data=Richness2)
leveneTest(aov(Richness~CO2*Temperature, data=Richness2))
# Checking for independence
durbinWatsonTest(model)

#><><><><><><><><><><><><><><><><><><><><><><><><>
#  Figure 2a -  Shared MOTUs across treatments
#><><><><><><><><><><><><><><><><><><><><><><><><>

source("SourceVenn.R")

### Shared OTUs
Share<-Richness2
Share<-as.data.table(Share)
Share2<-Share[,-c("Groups", "HeaderTank","Temperature","CO2","ARMS_Unit", "Richness")]
shareTreat<-sort(unique(Share2$Treatment))

share_list = list()
All_list<-list()

# Due to uneven number of ARMS per treatment taking a random draw of 5 units
  for (i in 1:length(shareTreat)){
    Treat<-subset(Share2, Treatment==shareTreat[i])
    Treat2<-sample_n(Treat, 5)
    Treat3<-melt(Treat2, id =c("Treatment", "Sample"), variable.name = "ID") 
    Treat4<-subset(Treat3, value > 0)
    Treat5<-ddply(Treat4, .(ID, Treatment), summarize,"Sum"=sum(value) )
    Treat6<-Treat5[,c(1)]
    share_list[[i]] <- sort(Treat6)
  }
  fill<- c("Control","Acidified","Heated","AcidifiedHeated")
  names(share_list) <- fill
  OLlist <- overLapper(setlist=share_list, sep="_", type="vennsets")
  counts <- sapply(OLlist$Venn_List, length)
  
  #Venn Bar Plot
  olBarplot(OLlist=OLlist, horiz=T, las=1, cex.names=0.6, main="Venn Bar Plot") 
  # brought into illustrator to alter colorations given the default coloration 
  # in the olBarplot function being based on OLlist$Complexity_Levels
 
  
#><><><><><><><><><><><><><><><><><><><><><><><><>
#  Figure 2b -  Treatment Boxplots
#><><><><><><><><><><><><><><><><><><><><><><><><>

TreatOrder<-c("Control", "Acidified","Heated","AcidifiedHeated")
Richness2$Treatment<-factor(Richness2$Treatment, levels = TreatOrder)

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

cols <- c("Control" = "#56B4E9", "Acidified" = "#F0E442",
          "Heated" = "#E69F00","AcidifiedHeated" = "#999999")

ggplot(Richness2, aes(x=Treatment, y=Richness)) + 
  geom_boxplot()+
  geom_point(aes(color = Treatment, fill = Treatment), pch = 21, size = 6, col = "Black", lwd=3)+
  scale_fill_manual(values=cols,aesthetics = c("Treatment", "fill"))+
  ylab("Observed MOTU Richness")+
  theme_bw()+
  xlab("")+
  ylim(0,130)+
  stat_summary(fun.data=data_summary)+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text=element_text(size=12))


#><><><><><><><><><><><><><><><><><><><><><><><><>
#  Beta Diversity
#><><><><><><><><><><><><><><><><><><><><><><><><>

Beta<-Richness2
Beta$Richness<-NULL

Beta2<-Beta[,8:ncol(Beta)]
metaBeta<-Beta[,1:7]


#><><><><><><><><><><><><><><>
#  PERMANOVA - Table S9
#><><><><><><><><><><><><><><>

# Community Composition - Presence/Absence

#first making dataframe binary
BetaJ<-decostand(Beta2, method = "pa")
# creating Jaccard dissimilarity matrix
BetaJ2<-vegdist(BetaJ, method = "jaccard")
adonis(BetaJ2 ~ Temperature*CO2,  data = metaBeta, permutations = 9999)
#Treatment groups for pairwise comparisons and PERMDISP 
groupTreat<-as.factor(Beta$Treatment)
# Treatment pairwise comparisons 
adonis.pair(BetaJ2,groupTreat, nper = 1000, corr.method = "fdr") #Table S10


# Community Structure - Relative Abudance 

# Hellinger transformation - convert to relative abundance and take sqrt
BetaRel<-decostand(Beta2, method = "hellinger")
BetaRel2<-vegdist(BetaRel, method = "bray")
adonis(BetaRel2 ~ Temperature*CO2,data = metaBeta,  permutations = 9999)
# Treatment pairwise comparisons
adonis.pair(BetaRel2,groupTreat, nper = 1000, corr.method = "fdr") #Table S10


#><><><><><><><><><><><<><><>
#  Figure 2C - PCoA - Jaccard 
#><><><><><><><><<><><><<><><>
Jac<-cmdscale(vegdist(BetaJ, method = "jaccard"), k = 2, eig = T, add = T)
labsJ<-paste("PCoA",1:4,"(", round(100*Jac$eig/sum(Jac$eig),2),"%)")
Jac1<-ordiplot(Jac, display="sites", xlab =labsJ[1],ylab=labsJ[2])
points(Jac1, "sites", pch = 19, col = "#56B4E9", select = Beta$Treatment == "Control")
points(Jac1, "sites", pch = 19, col = "#F0E442", select = Beta$Treatment == "Acidified")
points(Jac1, "sites", pch = 19, col = "#E69F00", select = Beta$Treatment == "Heated")
points(Jac1, "sites", pch = 19, col = "#999999", select = Beta$Treatment == "AcidifiedHeated")
ordispider(Jac1, Beta$Treatment, col=cols)    # col= usage needs vegan 2.4-0
ordiellipse(Jac1, Beta$Treatment, col=cols, draw="poly", conf=0.95,kind="se", label = F)


#><><><><><><><><><><><><><><<><><><<><><>
#  Figure 2D - PCoA - Bray-Curtis
#><><><><><><><><<><><><><><><<><><><<><><>

Bray<-cmdscale(vegdist(BetaRel, method = "bray"), k = 2, eig = T, add = T)
labsB<-paste("PCoA",1:4,"(", round(100*Bray$eig/sum(Bray$eig),2),"%)")
Bray1<-ordiplot(Bray, display="sites", xlab =labsB[1],ylab=labsB[2])
points(Bray1, "sites", pch = 19, col = "#56B4E9", select = Beta$Treatment == "Control")
points(Bray1, "sites", pch = 19, col = "#F0E442", select = Beta$Treatment == "Acidified")
points(Bray1, "sites", pch = 19, col = "#E69F00", select = Beta$Treatment == "Heated")
points(Bray1, "sites", pch = 19, col = "#999999", select = Beta$Treatment == "AcidifiedHeated")
ordispider(Bray1, Beta$Treatment, col=cols)    # col= usage needs vegan 2.4-0
ordiellipse(Bray1, Beta$Treatment, col=cols, draw="poly", conf=0.95,kind="se", label = F)


#><><><><><><><><><><>
#  Figure S4 - PCoA
#><><><><><><><><<><><>
# temperature
Jac_Temp<-ordiplot(Jac, display="sites", xlab =labsJ[1],ylab=labsJ[2])
points(Jac_Temp, "sites", pch = 19, col = "black", select = Beta$Temperature == "YesCO2")
points(Jac_Temp, "sites", pch = 19, col = "tan", select = Beta$Temperature == "NoCO2")
ordispider(Jac_Temp, Beta$Temperature, col=c("tan","black")) 
ordiellipse(Jac_Temp, Beta$Temperature, col=c("tan","black"), draw="poly", conf=0.95,kind="se")

# CO2 
Jac_Acid<-ordiplot(Jac, display="sites",xlab =labsJ[1],ylab=labsJ[2])
points(Jac_Acid, "sites", pch = 19, col = "black", select = Beta$CO2 == "YesCO2")
points(Jac_Acid, "sites", pch = 19, col = "tan", select = Beta$CO2 == "NoCO2")
ordispider(Jac_Acid, Beta$CO2, col=c("tan","black")) 
ordiellipse(Jac_Acid, Beta$CO2, col=c("tan","black"), draw="poly", conf=0.95,kind="se")

# Temperature 
Bray_Temp<-ordiplot(Jac, display="sites", xlab =labsB[1],ylab=labsB[2])
points(Bray_Temp, "sites", pch = 19, col = "black", select = Beta$Temperature == "YesCO2")
points(Bray_Temp, "sites", pch = 19, col = "tan", select = Beta$Temperature == "NoCO2")
ordispider(Bray_Temp, Beta$Temperature, col=c("tan","black")) 
ordiellipse(Bray_Temp, Beta$Temperature, col=c("tan","black"), draw="poly", conf=0.95,kind="se")

# CO2 
Bray_Acid<-ordiplot(Jac, display="sites", xlab =labsB[1],ylab=labsB[2])
points(Bray_Acid, "sites", pch = 19, col = "black", select = Beta$CO2 == "YesCO2")
points(Bray_Acid, "sites", pch = 19, col = "tan", select = Beta$CO2 == "NoCO2")
ordispider(Bray_Acid, Beta$CO2, col=c("tan","black")) 
ordiellipse(Bray_Acid, Beta$CO2, col=c("tan","black"), draw="poly", conf=0.95,kind="se")


#><><><><><><><><><><><><><><><><><><><><><>
#  PERDISP - Table S9 and S10 and Figure S5
#><><><><><><><><><><><><><><><><><><><><><>

groupTemp<-as.factor(Beta$Temperature)
groupAcid<-as.factor(Beta$CO2)

# Community Composition - Presence/Absence
#Treatment
B_jac<-betadisper(BetaJ2, groupTreat, bias.adjust = T)
permutest(B_jac)
permutest(B_jac, pairwise = TRUE)
TukeyHSD(B_jac)
boxplot(B_jac, col = c("light blue","yellow2","orange","gray"), ylim=c(0.20,0.80))
# Temperature
B_jacTemp<-betadisper(BetaJ2, groupTemp, bias.adjust = T)
anova(B_jacTemp)
boxplot(B_jacTemp,col = c("tan","gray22"), ylim=c(0.20,0.60))
# CO2
B_jacAcid<-betadisper(BetaJ2, groupAcid, bias.adjust = T)
anova(B_jacAcid)
boxplot(B_jacAcid,col = c("tan","gray22"), ylim=c(0.20,0.60))

#Double checking PCoA graph results
labs<-paste("PCoA",1:4,"(", round(100*B_jac$eig/sum(B_jac$eig),2),"%)")
plot(B_jac, main = "", hull = F, ellipse = T, 
     xlab =labs[1],ylab=labs[2], col = cols)

# Community Structure - Relative Abudance 

#Treatment
B_bray<-betadisper(BetaRel2, groupTreat, bias.adjust = T)
permutest(B_bray)
permutest(B_bray, pairwise = T)
TukeyHSD(B_bray)
boxplot(B_bray,col = c("light blue","yellow2","orange","gray"), ylim=c(0.20,0.80))
plot(TukeyHSD(B_bray))
# Temperature
B_brayTemp<-betadisper(BetaRel2, groupTemp, bias.adjust = T)
permutest(B_brayTemp)
boxplot(B_brayTemp,col = c("tan","gray22"), ylim=c(0.20,0.60))
# CO2
B_brayAcid<-betadisper(BetaRel2, groupAcid, bias.adjust = T)
permutest(B_brayAcid)
boxplot(B_brayAcid, col = c("tan","gray22"), ylim=c(0.20,0.60))


